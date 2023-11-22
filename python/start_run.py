import os
import shutil
import sys
from argparse import ArgumentParser
from datetime import datetime
from re import search, sub, split
from subprocess import run

import h5py
import numpy as np
from numpy import pi, zeros, append, sort, newaxis, digitize, stack, histogram, identity, diff, arange, \
	empty, maximum, minimum, cumsum, float64, average
from pandas import Series, Timestamp, notnull

from data_io import load_pulse_shape, parse_gas_components, load_beam_profile, load_inputs_table, \
	load_outputs_table, log_message, fill_in_template, write_row_to_outputs_table
from material import Material, get_solid_material_from_name, get_gas_material_from_components, isotope_symbol
from utilities import degrade_laser_pulse, gradient, select_key_indices, rebin


def start_run(code: str, name: str, stopping_power_mode: int, force: bool) -> None:
	""" arrange the inputs for a LILAC or IRIS run in a machine-readable format and submit it to slurm
	    :param code: one of "LILAC" or "IRIS"
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	    :param stopping_power_mode: one of 0 for no stopping, 1 for Li-Petrasso-Zylstra, or 2 for Maynard-Deutsch
	    :param force: whether to run this even if another version of this run already exists TODO: let the user hit "y" so they don't have to redo the command
	"""
	if name.endswith("*"):
		found_any_matches = False
		for full_name in load_inputs_table().index:
			if full_name[:len(name) - 1] == name[:-1]:
				print(f"starting {code} run for {full_name}...")
				start_run(code, full_name, stopping_power_mode, force)
				found_any_matches = True
		if not found_any_matches:
			raise ValueError(f"no rows in the run_inputs.csv table match '{name}'")

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
			print(f"this run seems to already be {current_status}.  to overwrite it, use the --force "
			      f"command line option.")
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
	    :raises KeyError: if it can't find the relevant information in the input table
	"""
	directory = f"runs/{name}/{code.lower()}"
	# get the run inputs from the run input table
	inputs_table = load_inputs_table()
	try:
		inputs = inputs_table.loc[name]
	except KeyError:
		raise KeyError(f"please add shot '{name}' to the run_inputs.csv file.")
	if inputs.ndim != 1:
		raise KeyError(f"run '{name}' appears more than once in the run_inputs.csv file!")

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
	shell_layer_materials = [get_solid_material_from_name(name) for name in split(r"\s*\+\s*", inputs["shell material"])]
	shell_layer_thicknesses = [float(x) for x in split(r"\s*\+\s*", inputs["shell thickness"])]

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
		fill_material: Material, shell_layer_materials: list[Material],
		shell_layer_thicknesses: list[float]) -> str:
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
			"fill protium percentage": fill_material.protium_fraction*100,
			"fill tritium percentage": fill_material.tritium_fraction*100,
			"fill EOS option": fill_material.eos,
			"fill opacity option": fill_material.opacity,
			"fill ionization option": fill_material.ionization,
			"fill pressure": fill_material.pressure,
			"fill radius": inputs['outer diameter']/2 - sum(shell_layer_thicknesses),
			# second prof namelist (shell)
			"shell thickness": shell_layer_thicknesses,
			"shell num cells": cell_counts,
			"shell material code": [material.material_code for material in shell_layer_materials],
			"shell protium percentage": [material.protium_fraction*100 for material in shell_layer_materials],
			"shell tritium percentage": [material.tritium_fraction*100 for material in shell_layer_materials],
			"shell EOS option": [material.eos for material in shell_layer_materials],
			"shell opacity option": [material.opacity for material in shell_layer_materials],
			"shell ionization option": [material.ionization for material in shell_layer_materials],
			"shell density": [material.density for material in shell_layer_materials],
			"shell density multiplier": inputs['shell density multiplier'] if notnull(inputs["shell density multiplier"]) else 1.,
			# third prof namelist (aluminum)
			"aluminum thickness": inputs['aluminum thickness'],
		},
		flags={
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
	"""
	# load information from the LILAC output file
	with h5py.File(f"runs/{name}/lilac/output.h5") as file:
		num_layers = file["target/material_name"].size
		material_codes = file["target/material_codes"][:, 0]
		material_names = file["target/material_name"][:]
		mean_atomic_masses = file["target/atomic_weight"][:]
		mean_atomic_numbers = average(file["target/component_nuclear_charge"][:, :],
		                              weights=file["target/component_abundance"], axis=1)

		# converting the component info to these fractions is a bit tricky
		nuclide_fractions = {nuclide: zeros(num_layers) for nuclide in ["H", "D", "T", "³He", "¹²C"]}
		for i in range(num_layers):
			for j in range(file["target/component_abundance"].shape[1]):
				nuclide = isotope_symbol(round(file["target/component_nuclear_charge"][i, j]),
				                         round(file["target/component_atomic_weight"][i, j]))
				if nuclide in nuclide_fractions:
					nuclide_fractions[nuclide][i] += file["target/component_abundance"][i, j]

		# the trickiest thing here is that we need to coarsen the spacio-temporal stuff
		burn_rate = gradient(file["time/total_dd_neutrons"][:] + file["time/total_dt_neutrons"][:],
		                     file["time/absolute_time"][:])
		time_bin_indices, key_time_indices = select_key_indices(burn_rate, 21)
		num_fine_times = burn_rate.size
		num_times = key_time_indices.size
		zone_height = diff(file["node/boundary_position"][:, :], axis=0)
		zone_mass = file["zone/mass_density"][:, 0]*zone_height[:, 0]
		key_node_indices, _ = select_key_indices(zone_mass, 50)
		# make sure layers don't get mixed together at the coarsened zone boundaries
		for i in range(num_layers):
			layer_boundary = file["target/last_zone"][i] - 1
			if layer_boundary not in key_node_indices:
				key_node_indices = append(key_node_indices, layer_boundary)
		key_node_indices = sort(key_node_indices)  # make sure it stays sorted
		num_fine_zones = zone_mass.shape[0]
		num_zones = key_node_indices.size - 1
		layer_index = digitize(key_node_indices[:-1], file["target/last_zone"][:] - 1)

		# then rebin the quantities to the new lower resolution
		time = file["time/absolute_time"][key_time_indices]
		node_position = file["node/boundary_position"][:, key_time_indices][key_node_indices, :]
		node_velocity = file["node/velocity"][:, key_time_indices][key_node_indices, :]
		zone_velocity = (node_velocity[0:-1] + node_velocity[1:])/2  # convert node velocity to zone velocity
		mass_density = rebin(file["zone/mass_density"][:, key_time_indices], key_node_indices,
		                     axis=0, weights=zone_height[:, key_time_indices])
		ion_temperature = rebin(file["zone/ion_temperature"][:, key_time_indices], key_node_indices,
		                        axis=0, weights=zone_height[:, key_time_indices])
		electron_temperature = rebin(file["zone/electron_temperature"][:, key_time_indices], key_node_indices,
		                             axis=0, weights=zone_height[:, key_time_indices])
		average_ionization = rebin(file["zone/average_z"][:, key_time_indices], key_node_indices,
		                           axis=0, weights=zone_height[:, key_time_indices])
		average_ionization2 = rebin(file["zone/average_z2"][:, key_time_indices], key_node_indices,
		                            axis=0, weights=zone_height[:, key_time_indices])
		material_fraction = identity(num_layers)[:, layer_index]  # this should be a sort of identity matrix
		# reactions in particular are tricky
		reactions = {}
		for key, cumulative_reactions in [("DT", file["zone/cumulative_dtn_yield"][:, :]), ("DD", file["zone/cumulative_ddn_yield"][:, :])]:
			# first convert to actual cumulative reactions at each time because I find the indexing easier that way
			cumulative_reactions = cumsum(cumulative_reactions.astype(float64), axis=1)
			# reduce the cumulative reactions to the key time bins
			cumulative_reactions = (
				cumulative_reactions[:, maximum(time_bin_indices - 1, 0)] +
				cumulative_reactions[:, minimum(time_bin_indices, num_fine_times - 1)])/2
			# differentiate in time to get the number of reactions in each section
			reactions_per_zone = diff(cumulative_reactions, axis=1)
			# then rebin in space to the new coarser space bins
			reactions_per_bin = empty((num_zones, num_times))
			for j in range(num_times):
				reactions_per_bin[:, j] = histogram(
					arange(num_fine_zones), weights=reactions_per_zone[:, j], bins=key_node_indices)[0]
			reactions[key] = reactions_per_bin

	# save everything to the new directory
	directory = f"runs/{name}/iris-model{stopping_power_mode}"
	os.makedirs(f"{directory}/profiles", exist_ok=True)

	# write all of the new HDF5 input files
	keV = 1e3*1.602e-19
	for j in range(num_times):  # the filenames index from 1 because Fortran does
		with h5py.File(f"{directory}/profiles/time{j + 1}.h5", "w") as file:
			file["title"] = name
			file["n_r_cells"] = num_zones
			file["n_theta_cells"] = 1  # one theta cell because LILAC is 1D
			file["n_phi_cells"] = 1  # one phi cell because LILAC is 1D
			file["n_materials"] = num_layers
			file["time"] = time[j]  # ns
			file["material_names"] = material_names
			file["r"] = node_position[newaxis, newaxis, :, j]*1e-6  # m
			file["t"] = [0, pi]  # rad
			file["p"] = [0, 2*pi]  # rad
			file["mass_density"] = mass_density[newaxis, newaxis, :, j]/1e-3  # kg/m^3
			file["ion_temperature"] = ion_temperature[newaxis, newaxis, :, j]*keV  # J
			file["electron_temperature"] = electron_temperature[newaxis, newaxis, :, j]*keV  # J
			file["z_effective"] = (average_ionization2/average_ionization)[newaxis, newaxis, :, j]
			file["r_velocity"] = zone_velocity[newaxis, newaxis, :, j]*1e-2  # m/s
			file["t_velocity"] = zeros((1, 1, num_zones))
			file["p_velocity"] = zeros((1, 1, num_zones))
			file["material_fraction"] = material_fraction[:, newaxis, newaxis, :]
			file["reactions"] = stack([reactions["DT"][newaxis, newaxis, :, j],
			                           reactions["DD"][newaxis, newaxis, :, j],
			                           zeros((1, 1, num_zones))], axis=0)

	# compose the input deck and bash script
	current_working_directory = os.getcwd().replace('\\', '/')
	input_deck = fill_in_template(
		"iris_input_deck.txt", {
			"submitted": datetime.today().isoformat(" "),
			"sanitized name": sub(r"[/\\{}<> ~#%&*?]", "-", name),
			"number of times": time.size,
			"fill material code": material_codes[0],
			"fill hydrogen fraction": nuclide_fractions['H'][0],
			"fill deuterium fraction": nuclide_fractions['D'][0],
			"fill tritium fraction": nuclide_fractions['T'][0],
			"fill helium3 fraction": nuclide_fractions['³He'][0],
			"fill carbon fraction": nuclide_fractions['¹²C'][0],
			"fill mean atomic mass": mean_atomic_masses[0]*1.66054e-27,  # kg
			"fill max ionization": mean_atomic_numbers[0],
			"shell material code": material_codes[1],
			"shell hydrogen fraction": nuclide_fractions['H'][1],
			"shell deuterium fraction": nuclide_fractions['D'][1],
			"shell tritium fraction": nuclide_fractions['T'][1],
			"shell helium3 fraction": nuclide_fractions['³He'][1],
			"shell carbon fraction": nuclide_fractions['¹²C'][1],
			"shell mean atomic mass": mean_atomic_masses[1]*1.66054e-27,  # kg
			"shell max ionization": mean_atomic_numbers[1],
			"stopping power model": stopping_power_mode,
			"directory": f"{current_working_directory}/{directory}",
		},
		loops = {"j": [str(j + 1) for j in range(num_times)]},  # index from 1 because IRIS does
	)
	with open(f"{directory}/inputdeck.txt", "w") as file:
		file.write(input_deck)

	bash_script = fill_in_template(
		"iris_run.sh", {
			"name": name,
			"root": current_working_directory,
			"folder": f"iris-model{stopping_power_mode}"
		})
	with open(f"{directory}/run.sh", "w") as file:
		file.write(bash_script),
	return f"{directory}/run.sh"


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
		help="the number of the plasma stopping-power model to use for charged particles (0 = none, 1 = Li-Petrasso-Zylstra, 2 = Maynard-Deutsch)",)
	parser.add_argument(
		"--force", action="store_true",
		help="whether to overwrite any previous iterations of this run")
	args = parser.parse_args(sys.argv[2:])  # TODO: add arguments to override flux limiter, density, and laser degradation

	try:
		start_run(code, args.name, args.stopping_mode, args.force)
		sys.exit(0)
	except Exception as e:
		print(f"Error!", *e.args)
		sys.exit(1)
