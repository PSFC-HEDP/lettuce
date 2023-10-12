import os
import shutil
from argparse import ArgumentParser
from datetime import datetime
from re import sub, DOTALL
from subprocess import call

import numpy as np
from pandas import Series, Timestamp, concat

from python.utilities import load_pulse_shape, get_shell_material_from_name, get_gas_material_from_components, \
	parse_gas_components, load_beam_profile, Material, load_inputs_table, load_outputs_table, log_message, \
	submit_slurm_job


def start_lilac_run(name: str, force: bool) -> None:
	""" arrange the inputs for a LILAC run in a machine-readable format and submit it to slurm
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	    :param force: whether to run this
	"""
	# start by loading up the run input and output tables
	inputs_table = load_inputs_table()
	outputs_table = load_outputs_table()

	# assess the current state of this run
	if (name, "LILAC") not in outputs_table.index:
		current_status = "new"
	elif not os.path.isdir(f"runs/{name}/lilac"):
		current_status = "new"
	else:
		current_status = outputs_table.loc[(name, "LILAC"), "status"]
	# if it's already running and the user didn't force it, stop work immediately
	if current_status in ["running", "completed", "pending"] and not force:
		print(f"this run seems to already be {current_status}.  to overwrite it, use the --force "
		      f"command line option.")
		return
	# if it's already running and the user did force it, cancel the run
	elif current_status == "pending" or current_status == "running":
		print(f"cancelling the current run...")
		call(["scancel", outputs_table.loc[(name, "LILAC"), "slurm ID"]])
	# if it's already run and the user did force it, clear the previous output
	elif current_status == "completed":
		print(f"overwriting the current run...")
		shutil.rmtree(f"runs/{name}/lilac")

	# get the run inputs from the run input table
	try:
		inputs = inputs_table.loc[name]
	except KeyError:
		print(f"please add shot '{name}' to the run_inputs.csv file.")
		return
	if inputs.ndim != 1:
		print(f"run '{name}' appears more than once in the run_inputs.csv file!")
		return

	# load the laser data
	try:
		pulse_time, pulse_power = load_pulse_shape(inputs["pulse shape"], inputs["laser energy"])
		beam_radius, beam_intensity = load_beam_profile(inputs["beam profile"])
	except (IOError, ValueError) as e:
		print(e)
		return

	# parse the materials
	shell_material = get_shell_material_from_name(inputs["shell material"])
	fill_material = get_gas_material_from_components(parse_gas_components(inputs["fill"]))

	# compose the input deck and bash script
	input_deck = build_lilac_input_deck(name, inputs, fill_material, shell_material)
	bash_script = build_lilac_bash_script(name)

	# save all of the inputs and the script to the run directory
	os.makedirs(f"runs/{name}/lilac", exist_ok=True)
	with open(f"runs/{name}/lilac/input_deck.txt", "w") as file:
		file.write(input_deck)
	with open(f"runs/{name}/lilac/run_lilac.sh", "w") as file:
		file.write(bash_script)
	np.savetxt(f"runs/{name}/lilac/pulse_shape.txt",
	           np.stack([pulse_time, pulse_power], axis=1), delimiter=",")
	np.savetxt(f"runs/{name}/lilac/beam_profile.txt",
	           np.stack([beam_radius, beam_intensity], axis=1), delimiter=" ")

	# finally, submit the slurm job
	slurm_ID = submit_slurm_job("run_lilac.sh")

	# update our records
	if (name, "LILAC") not in outputs_table.index:
		outputs = Series(
			name=(name, "LILAC"), data={
				"status":          "pending",
				"status changed":  Timestamp.now(),
				"slurm ID":        slurm_ID})
		outputs_table = concat([outputs_table, outputs.to_frame().T])
	outputs_table.to_csv("run_outputs.csv")
	log_message(f"LILAC run '{name}' (slurm ID {slurm_ID}) is submitted to slurm.")


def build_lilac_input_deck(
		name: str, inputs: Series, fill_material: Material, shell_material: Material) -> str:
	""" construct the input deck corresponding to the given inputs and return it as a str """
	with open("resources/templates/lilac_input_deck.txt") as template_file:
		input_deck = template_file.read()

	# rhydro namelist
	input_deck = sub("<<name>>", name, input_deck)
	input_deck = sub("<<absorption fraction>>", f"{inputs['absorption fraction']:.5f}", input_deck)
	input_deck = sub("<<nonthermal model>>", "none" if inputs["flux limiter"] != 0 else "", input_deck)
	input_deck = sub("<<flux limiter>>", f"{inputs['flux limiter']:.3f}", input_deck)

	# first prof namelist (fill)
	input_deck = sub("<<fill material code>>", f"{fill_material.material_code:d}", input_deck)
	input_deck = sub("<<fill protium fraction>>", f"{fill_material.protium_fraction:.5f}", input_deck)
	input_deck = sub("<<fill tritium fraction>>", f"{fill_material.tritium_fraction:.5f}", input_deck)
	input_deck = sub("<<fill pressure>>", f"{fill_material.pressure:.4f}", input_deck)
	input_deck = sub(
		"<<fill radius>>", f"{inputs['outer diameter']/2 - inputs['shell thickness']:.3f}", input_deck)

	# second prof namelist (shell)
	input_deck = sub("<<shell material code>>", f"{shell_material.material_code:d}", input_deck)
	input_deck = sub("<<shell protium fraction>>", f"{shell_material.protium_fraction:.5f}", input_deck)
	input_deck = sub("<<shell tritium fraction>>", f"{shell_material.tritium_fraction:.5f}", input_deck)
	input_deck = sub("<<shell thickness>>", f"{inputs['shell thickness']:.3f}", input_deck)
	if shell_material.density is not None:
		input_deck = sub("<<(end)?if density specified>>\n", "", input_deck)
		input_deck = sub("<<shell density>>", f"{inputs['shell density']:.4f}", input_deck)
	else:
		input_deck = sub("<<if density specified>>.*<<endif density specified>>\n", "", input_deck, flags=DOTALL)
	input_deck = sub("<<shell density multiplier>>", f"{inputs['density multiplier']:.5f}", input_deck)

	# third prof namelist (aluminum)
	if inputs["aluminum thickness"] > 0:
		input_deck = sub("<<(end)?if aluminum>>\n", "", input_deck)
		input_deck = sub("<<aluminum thickness>>", f"{inputs['aluminum thickness']:.3f}", input_deck)
	else:
		input_deck = sub("<<if aluminum>>.*<<endif aluminum>>\n", "", input_deck, DOTALL)

	input_deck = sub("<<submitted>>", datetime.today().isoformat(" "), input_deck)
	return input_deck


def build_lilac_bash_script(name: str) -> str:
	""" construct the bash script that will run the given lilac and return it as a str """
	with open("resources/templates/run_lilac.sh") as template_file:
		bash_script = template_file.read()

	base_directory = os.getcwd().replace('\\', '/')
	bash_script = sub("<<name>>", name, bash_script)
	bash_script = sub("<<directory>>", f"{base_directory}/runs/{name}/lilac", bash_script)
	return bash_script


if __name__ == "__main__":
	parser = ArgumentParser(
		prog="start_lilac_run.sh",
		description = "Arrange the inputs for a LILAC run in a machine-readable format and submit it to slurm")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	parser.add_argument(
		"--force", action="store_true",
		help="whether to overwrite any previous iterations of this run")
	args = parser.parse_args()

	start_lilac_run(args.name, args.force)
