import os
from os import path
import shutil
import sys
from argparse import ArgumentParser
from datetime import datetime
from re import search, sub, split
from subprocess import run
from typing import List, Optional

import numpy as np
from pandas import Series, Timestamp, notnull

from data_io import load_pulse_shape, parse_gas_components, load_beam_profile, load_inputs_table, \
	load_outputs_table, log_message, fill_in_template, write_row_to_outputs_table
from material import Material, get_solid_material_from_name, get_gas_material_from_components, \
	parse_isotope_symbol
from utilities import InvalidSimulationError, RecordNotFoundError, degrade_laser_pulse


def start_run(code: str, name: str, stopping_power_option: Optional[int], force: bool) -> None:
	""" arrange the inputs for a LILAC or IRIS run in a machine-readable format and submit it to slurm
	    :param code: one of "LILAC" or "IRIS"
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	    :param stopping_power_option: one of 0 for no stopping, 1 for Li-Petrasso-Zylstra, or 2 for Maynard-Deutsch.
	                                  if none is specified, it will use 1 but won't specify that in the filenames.
	    :param force: whether to run this even if another version of this run already exists TODO: let the user hit "y" so they don't have to redo the command
	    :raise RecordNotFoundError: if it can't find the relevant simulation inputs to build the input deck
	    :raise InvalidSimulationError: if for whatever reason the relevant simulation inputs won't work
	"""
	# expand any star expressions in the name
	if name.endswith("*"):
		found_any_matches = False
		for full_name in load_inputs_table().index:
			if full_name[:len(name) - 1] == name[:-1]:
				print(f"starting {code} run for {full_name}...")
				start_run(code, full_name, stopping_power_option, force)
				found_any_matches = True
		if not found_any_matches:
			raise RecordNotFoundError(f"no rows in the run_inputs.csv table match '{name}'")
		else:
			return

	# augment the name if necessary to incorporate any command-line options
	# TODO: add command-line options for things like flux limiter and pulse degradation
	lilac_name = name
	if stopping_power_option is not None:
		if code == "LILAC":
			raise InvalidSimulationError("you can't set the stopping power model for LILAC; that's an IRIS option.")
		stopping_power_model_number = stopping_power_option
		name = lilac_name + f"/model{stopping_power_model_number:d}"
	else:
		stopping_power_model_number = 1
		name = lilac_name

	if code == "LILAC":
		# assess the current state of this run
		outputs_table = load_outputs_table()
		if name not in outputs_table.index:
			current_lilac_status = "new"
		elif not path.isdir(f"runs/{name}"):
			current_lilac_status = "new"
		else:
			current_status: str = outputs_table.loc[name, "status"]
			if current_status.endswith(" LILAC"):
				current_lilac_status = current_status.split()[0]
			elif current_status.endswith(" IRIS"):
				current_lilac_status = "completed"
			else:
				raise ValueError(f"I don't think '{current_status}' is a valid status.")

		# if it's already running and the user didn't force it, stop work immediately
		if current_lilac_status in ["running", "completed", "pending"] and not force: # TODO: it seems like IRIS could also benefit from this part of the code
			answer = input(
				f"This run seems to already be {current_lilac_status} (use `--force` to ignore this warning).  "
				f"Would you like to overwrite it?  [y/N] "
			).lower()
			if answer == "" or answer == "n" or answer == "no":
				print("Mm-hmm, I thought not.")
				return
			elif answer == "y" or answer == "yes":
				print("you got it, boss.")
			else:
				print(f"I'm not sure what '{answer}' is supposed to mean so I'm going to assume that's a 'no'.")
				return
		# if it's already running and the user did force it, cancel the run
		elif current_lilac_status == "pending" or current_lilac_status == "running":
			print(f"cancelling the current run...")
			try:
				run(["scancel", outputs_table.loc[name, "slurm ID"]])
			except FileNotFoundError:
				raise FileNotFoundError("I can't call scancel. are you sure slurm is installed and loaded?")
		# if it's already run and the user did force it, clear the previous output
		elif current_lilac_status == "completed":
			print(f"overwriting the previous run...")
			shutil.rmtree(f"runs/{name}")

	if code == "LILAC":
		script_path = prepare_lilac_inputs(name)
	elif code == "IRIS":
		script_path = prepare_iris_inputs(lilac_name, name, stopping_power_model_number)
	else:
		print(f"unrecognized code: '{code}'")
		return

	# finally, submit the slurm job
	try:
		submission = run(["sbatch", script_path],
		                 check=True, capture_output=True, text=True)
	except FileNotFoundError:
		raise FileNotFoundError("I can't call sbatch. are you sure slurm is installed and loaded?")
	slurm_ID = search(r"batch job ([0-9]+)", submission.stdout).group(1)

	# update our records
	write_row_to_outputs_table({
		"name": name,
		"status": f"pending {code}",
		"status changed": Timestamp.now(),
		"slurm ID": slurm_ID,
	})
	log_message(f"{code} run '{name}' is submitted to slurm (slurm ID {slurm_ID}).")


def prepare_lilac_inputs(name: str) -> str:
	""" load the LILAC input parameters from run_inputs.csv and save the pulse shape, beam profile,
	    and input deck to the relevant directory.
	    :raise RecordNotFoundError: if it can't find the relevant information in the input table
	"""
	# get the run inputs from the run input table
	inputs_table = load_inputs_table()
	try:
		inputs = inputs_table.loc[name]
	except KeyError:
		raise RecordNotFoundError(f"please add shot '{name}' to the run_inputs.csv file.")
	if inputs.ndim != 1:
		raise RecordNotFoundError(f"run '{name}' appears more than once in the run_inputs.csv file!")

	# save all of the inputs and the script to the run directory
	os.makedirs(f"runs/{name}", exist_ok=True)

	# load the laser data
	pulse_time, pulse_power = load_pulse_shape(inputs["pulse shape"], inputs["laser energy"])
	if notnull(inputs["laser degradation"]):
		pulse_power = degrade_laser_pulse(pulse_power, inputs["laser degradation"])
	np.savetxt(f"runs/{name}/pulse_shape.txt",
	           np.stack([pulse_time, pulse_power], axis=1), delimiter=",")  # type: ignore

	beam_radius, beam_intensity = load_beam_profile(inputs["beam profile"])
	np.savetxt(f"runs/{name}/beam_profile.txt",
	           np.stack([beam_radius, beam_intensity], axis=1),  # type: ignore
	           delimiter=" ", fmt="%.6f")  # type: ignore

	# parse the materials
	fill_material = get_gas_material_from_components(parse_gas_components(inputs["fill"]))
	shell_layer_materials = [get_solid_material_from_name(name) for name in split(r"\s*/\s*", inputs["shell material"])]
	shell_layer_thicknesses = [float(x) for x in split(r"\s*/\s*", inputs["shell thickness"])]

	pulse_end_time = pulse_time[pulse_power >= 1/2*np.max(pulse_power)][-1]

	# compose the input deck
	input_deck = build_lilac_input_deck(name, inputs, pulse_end_time, fill_material,
	                                    shell_layer_materials, shell_layer_thicknesses)
	with open(f"runs/{name}/lilac_data_input.txt", "w") as file:
		file.write(input_deck)

	# compose the bash script
	bash_script = build_lilac_bash_script(name)
	with open(f"runs/{name}/lilac.sh", "w") as file:
		file.write(bash_script)

	return f"runs/{name}/lilac.sh"


def build_lilac_input_deck(
		name: str, inputs: Series, laser_off_time: float,
		fill_material: Material, shell_layer_materials: List[Material],
		shell_layer_thicknesses: List[float]) -> str:
	""" construct the input deck corresponding to the given inputs and return it as a str """
	# the one thing we need to calculate ourselves is the number of cells in the shell
	cell_counts = []
	for i in range(len(shell_layer_materials)):
		if shell_layer_materials[i].density is not None:
			areal_density = shell_layer_materials[i].density*shell_layer_thicknesses[i]
		else:
			areal_density = 1.0*shell_layer_thicknesses[i]
		if notnull(inputs["shell density multiplier"]):
			areal_density *= inputs["shell density multiplier"]
		cell_counts.append(max(10, min(500, round(areal_density*15))))

	# we also need to enumerate the component species of the fuel if it's a custom material
	if fill_material.material_code < 0:
		custom_fill = True
		component_symbols = list(fill_material.components.keys())
		component_weights, component_charges = [], []
		component_codes, component_abundances = [], []
		for symbol in component_symbols:
			charge, weight = parse_isotope_symbol(symbol)
			component_weights.append(weight)
			component_charges.append(charge)
			component_codes.append(charge if symbol not in "DT" else 21)
			component_abundances.append(fill_material.components[symbol])
	else:
		custom_fill = False
		component_symbols, component_abundances = None, None
		component_weights, component_charges = None, None
		component_codes = None

	return fill_in_template(
		"lilac_input_deck.txt",
		parameters={
			# header
			"submitted": datetime.today().isoformat(" "),
			# rhydro namelist
			"sanitized name": sub(r"[/\\{}<> ~#%&*?]", "-", name),
			"absorption fraction": inputs['absorption fraction'],
			"nonthermal model": "none" if inputs["flux limiter"] > 0 else "vgon",
			"flux limiter": inputs['flux limiter'],
			"laser off time": laser_off_time,
			# first prof namelist (fill)
			"fill material code": fill_material.material_code,
			"fill protium percentage": fill_material.components.get("H", 0)*100,
			"fill tritium percentage": fill_material.components.get("T", 0)*100,
			"fill EOS option": fill_material.eos,
			"fill opacity option": fill_material.opacity,
			"fill ionization option": fill_material.ionization,
			"fill pressure": fill_material.pressure,
			"fill radius": inputs['outer diameter']/2 - sum(shell_layer_thicknesses),
			# mater namelist (fill mix)
			"true fill material code": abs(fill_material.material_code) if custom_fill else None,
			"fill opacity table": fill_material.opacity_table,
			"fill component symbols": component_symbols,
			"fill component abundances": component_abundances,
			"fill component atomic weights": component_weights,
			"fill component atomic numbers": component_charges,
			"fill component material codes": component_codes,
			# second prof namelist (shell)
			"shell thickness": shell_layer_thicknesses,
			"shell num cells": cell_counts,
			"shell material code": [material.material_code for material in shell_layer_materials],
			"shell protium percentage": [material.components.get("H", 0)*100 for material in shell_layer_materials],
			"shell tritium percentage": [material.components.get("T", 0)*100 for material in shell_layer_materials],
			"shell EOS option": [material.eos for material in shell_layer_materials],
			"shell opacity option": [material.opacity for material in shell_layer_materials],
			"shell ionization option": [material.ionization for material in shell_layer_materials],
			"shell density": [material.density for material in shell_layer_materials],
			"shell density multiplier": inputs['shell density multiplier'] if notnull(inputs["shell density multiplier"]) else 1.,
			# third prof namelist (aluminum)
			"aluminum thickness": inputs['aluminum thickness'],
		},
		flags={
			"custom fill": custom_fill,
			"density specified": [material.density is not None for material in shell_layer_materials],
			"aluminum": inputs["aluminum thickness"] > 0,
		},
		loops={
			"i": range(len(shell_layer_materials)),
		},
	)


def build_lilac_bash_script(name: str) -> str:
	""" construct the bash script that will run the given lilac and return it as a str """
	return fill_in_template(
		"lilac_run.sh", {
			"name": name,
			"basename": path.basename(name),
			"root": os.getcwd().replace('\\', '/'),
		})


def prepare_iris_inputs(lilac_name: str, iris_name: str, stopping_power_model_number: int) -> str:
	""" load the relevant LILAC outputs and save them and a matching input deck to the relevant directory.
	    :param lilac_name: the name of the LILAC run off of which we're basing this
	    :param iris_name: the name of the IRIS run to set up
	    :param stopping_power_model_number: one of 0 for no stopping, 1 for Li-Petrasso-Zylstra, or 2 for Maynard-Deutsch.
	    :return: the filepath of the prepared bash script
	    :raise RecordNotFoundError: if the LILAC input can't be found
	    :raise InvalidSimulationError: if the shot has no DT reactions
	"""
	import lotus

	current_working_directory = os.getcwd().replace('\\', '/')

	# use Lotus to generate the input deck and profile HDF5 files
	try:
		lilac_solution = lotus.lilac.LilacSolution(hdf5_file_path=f"runs/{lilac_name}/lilac_output_{path.basename(lilac_name)}.h5")
	except FileNotFoundError:
		raise RecordNotFoundError(f"There does not appear to be any LILAC output in runs/{iris_name}")
	if np.all(lilac_solution.burn.reaction_rate("D(T,n)He4") == 0):
		raise InvalidSimulationError(f"This simulation has no DT reactions.  IRIS only works on DT shots.")
	lotus.postprocessors.iris.IRISInputDeck(
		shot_number=lilac_solution.shot_number,
		hydrocode_solution=lilac_solution,
		output_file=f"{current_working_directory}/runs/{iris_name}/inputdeck.txt",
	)

	# make some adjustments (these can probably be done with Lotus but I don't have access to the documentation rite now so)
	with open(f"runs/{iris_name}/inputdeck.txt", "r") as file:
		input_deck = file.read()
	# add the charged particle transport model option
	input_deck = sub(
		r"&scatter([^/]*)/",
		f"&scatter\\1\n\n"
		f"    ! stopping power model (0 = no stopping, 1 = Li-Petrasso-Zylstra, 2 = Maynard-Deutsch)\n"
		f"    charged_particle_transport_model = {stopping_power_model_number:d}\n"
		f"/",
		input_deck
	)
	with open(f"runs/{iris_name}/inputdeck.txt", "w") as file:
		file.write(input_deck)

	bash_script = fill_in_template(
		"iris_run.sh", {
			"name": iris_name,
			"basename": path.basename(iris_name),
			"root": current_working_directory,
		})
	with open(f"runs/{iris_name}/iris.sh", "w") as file:
		file.write(bash_script),
	return f"runs/{iris_name}/iris.sh"


if __name__ == "__main__":
	code = sys.argv[1]
	parser = ArgumentParser(
		prog=f"start_{code.lower()}_run.sh",
		description=f"Arrange the inputs for a {code} run in a machine-readable format and submit it to slurm")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	parser.add_argument(
		"--stopping_model", type=int, default=None,
		help="the number of the plasma stopping-power model to use for charged particles "
		     "(0 = none, 1 = Li-Petrasso-Zylstra, 2 = Maynard-Deutsch)")
	parser.add_argument(
		"--force", action="store_true",
		help="whether to overwrite any previous iterations of this run (if --force is not set and you try to start a "
		     "run that already exists, it will ask you if you're sure)")
	args = parser.parse_args(sys.argv[2:])  # TODO: add arguments to override flux limiter, density, and laser degradation

	try:
		start_run(code, args.name, args.stopping_model, args.force)
	except (RecordNotFoundError, InvalidSimulationError) as e:
		print("Error!", e)
		sys.exit(1)
	else:
		sys.exit(0)
