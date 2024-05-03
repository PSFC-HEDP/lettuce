import os
import shutil
import sys
from argparse import ArgumentParser
from datetime import datetime
from re import search, sub, split
from subprocess import run
from typing import List

import numpy as np
from pandas import Series, Timestamp, notnull

from data_io import load_pulse_shape, parse_gas_components, load_beam_profile, load_inputs_table, \
	load_outputs_table, log_message, fill_in_template, write_row_to_outputs_table
from material import Material, get_solid_material_from_name, get_gas_material_from_components, isotope_symbol, \
	parse_isotope_symbol
from utilities import InvalidSimulationError, RecordNotFoundError, degrade_laser_pulse, gradient, select_key_indices, rebin


def start_run(code: str, name: str, stopping_power_mode: int, force: bool) -> None:
	""" arrange the inputs for a LILAC or IRIS run in a machine-readable format and submit it to slurm
	    :param code: one of "LILAC" or "IRIS"
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	    :param stopping_power_mode: one of 0 for no stopping, 1 for Li-Petrasso-Zylstra, or 2 for Maynard-Deutsch
	    :param force: whether to run this even if another version of this run already exists TODO: let the user hit "y" so they don't have to redo the command
		:raise RecordNotFoundError: if it can't find the relevant simulation inputs to build the input deck
		:raise InvalidSimulationError: if for whatever reason the relevant simulation inputs won't work
	"""
	if name.endswith("*"):
		found_any_matches = False
		for full_name in load_inputs_table().index:
			if full_name[:len(name) - 1] == name[:-1]:
				print(f"starting {code} run for {full_name}...")
				start_run(code, full_name, stopping_power_mode, force)
				found_any_matches = True
		if not found_any_matches:
			raise KeyError(f"no rows in the run_inputs.csv table match '{name}'")
		else:
			return

	if code == "LILAC":
		# assess the current state of this run TODO: the outputs table needs to distinguish the LILAC state from the IRIS state, and warn you if you try to do IRIS without LILAC first
		outputs_table = load_outputs_table()
		if name not in outputs_table.index:
			current_status = "new"
		elif not os.path.isdir(f"runs/{name}/lilac"):
			current_status = "new"
		else:
			current_status = outputs_table.loc[name, "status"]
		# if it's already running and the user didn't force it, stop work immediately
		if current_status in ["running", "completed", "pending"] and not force:
			answer = input(
				f"This run seems to already be {current_status}.  Would you like to overwrite it (use the --force "
				f"command line option to ignore this warning)?  [y/N]"
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
		elif current_status == "pending" or current_status == "running":
			print(f"cancelling the current run...")
			try:
				run(["scancel", outputs_table.loc[name, "slurm ID"]])
			except FileNotFoundError:
				raise FileNotFoundError("I can't call scancel. are you sure slurm is installed and loaded?")
		# if it's already run and the user did force it, clear the previous output
		elif current_status == "completed":
			print(f"overwriting the previous run...")
			shutil.rmtree(f"runs/{name}/lilac")

	if code == "LILAC":
		script_path = prepare_lilac_inputs(name)
	elif code == "IRIS":
		script_path = prepare_iris_inputs(name, stopping_power_mode)
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
		"status": "pending",
		"status changed": Timestamp.now(),
		"slurm ID": slurm_ID,
	})
	log_message(f"{code} run '{name}' is submitted to slurm (slurm ID {slurm_ID}).")


def prepare_lilac_inputs(name: str) -> str:
	""" load the LILAC input parameters from run_inputs.csv and save the pulse shape, beam profile,
	    and input deck to the relevant directory.
	    :raise RecordNotFoundError: if it can't find the relevant information in the input table
	"""
	directory = f"runs/{name}/{code.lower()}"
	# get the run inputs from the run input table
	inputs_table = load_inputs_table()
	try:
		inputs = inputs_table.loc[name]
	except KeyError:
		raise RecordNotFoundError(f"please add shot '{name}' to the run_inputs.csv file.")
	if inputs.ndim != 1:
		raise RecordNotFoundError(f"run '{name}' appears more than once in the run_inputs.csv file!")

	# save all of the inputs and the script to the run directory
	os.makedirs(directory, exist_ok=True)

	# load the laser data
	pulse_time, pulse_power = load_pulse_shape(inputs["pulse shape"], inputs["laser energy"])
	if notnull(inputs["laser degradation"]):
		pulse_power = degrade_laser_pulse(pulse_power, inputs["laser degradation"])
	np.savetxt(f"{directory}/pulse_shape.txt",
	           np.stack([pulse_time, pulse_power], axis=1), delimiter=",")  # type: ignore

	beam_radius, beam_intensity = load_beam_profile(inputs["beam profile"])
	np.savetxt(f"{directory}/beam_profile.txt",
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
	with open(f"{directory}/lilac_data_input.txt", "w") as file:
		file.write(input_deck)

	# compose the bash script
	bash_script = build_lilac_bash_script(name)
	with open(f"{directory}/run.sh", "w") as file:
		file.write(bash_script)

	return f"{directory}/run.sh"


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
			"root": os.getcwd().replace('\\', '/'),
		})


def prepare_iris_inputs(name: str, stopping_power_mode: int) -> str:
	""" load the relevant LILAC outputs and save them and a matching input deck to the relevant directory.
	    :raise RecordNotFoundError: if the LILAC input can't be found
		:raise InvalidSimulationError: if the shot has no DT reactions
	"""
	import lotus

	lilac_directory = f"runs/{name}/lilac/"
	iris_directory = f"runs/{name}/iris-model{stopping_power_mode}" # TODO: I don't really need to separate the lilac and iris directories like this.
	current_working_directory = os.getcwd().replace('\\', '/')

	try:
		lilac_solution = lotus.lilac.LilacSolution(hdf5_file_path=f"{lilac_directory}/output.h5")
	except FileNotFoundError:
		raise RecordNotFoundError(f"There does not appear to be any LILAC output in {lilac_directory}")
	if np.all(lilac_solution.burn.reaction_rate("D(T,n)He4") == 0):
		raise InvalidSimulationError(f"This simulation has no DT reactions.  IRIS only works on DT shots.")
	lotus.postprocessors.iris.IRISInputDeck(
		shot_number=lilac_solution.shot_number,
		hydrocode_solution=lilac_solution,
		output_file=f"{iris_directory}/inputdeck.txt",
	)

	bash_script = fill_in_template(
		"iris_run.sh", {
			"name": name,
			"root": current_working_directory,
			"folder": f"iris-model{stopping_power_mode}"
		})
	with open(f"{iris_directory}/run.sh", "w") as file:
		file.write(bash_script),
	return f"{iris_directory}/run.sh"


if __name__ == "__main__":
	code = sys.argv[1]
	parser = ArgumentParser(
		prog=f"start_{code.lower()}_run.sh",
		description=f"Arrange the inputs for a {code} run in a machine-readable format and submit it to slurm")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	parser.add_argument(
		"--stopping_mode", type=int, default=1,
		help="the number of the plasma stopping-power model to use for charged particles "
		     "(0 = none, 1 = Li-Petrasso-Zylstra, 2 = Maynard-Deutsch)")
	parser.add_argument(
		"--force", action="store_true",
		help="whether to overwrite any previous iterations of this run (if --force is not set and you try to start a "
		     "run that already exists, it will ask you if you're sure)")
	args = parser.parse_args(sys.argv[2:])  # TODO: add arguments to override flux limiter, density, and laser degradation

	try:
		start_run(code, args.name, args.stopping_mode, args.force)
	except (RecordNotFoundError, InvalidSimulationError) as e:
		print("Error!", e)
		sys.exit(1)
	else:
		sys.exit(0)
