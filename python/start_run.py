import os
import shutil
import sys
from argparse import ArgumentParser
from datetime import datetime
from re import search, sub
from subprocess import run

import h5py
import numpy as np
from numpy import pi, zeros, append, sort, newaxis, digitize, stack, histogram, identity, diff, arange, \
	empty, maximum, minimum, cumsum, float64, average
from pandas import Series, Timestamp, notnull

from data_io import load_pulse_shape, parse_gas_components, load_beam_profile, load_inputs_table, \
	load_outputs_table, log_message, fill_in_template, write_row_to_outputs_table
from material import Material, get_solid_material_from_name, get_gas_material_from_components, nuclide_symbol
from utilities import degrade_laser_pulse, gradient, select_key_indices, rebin


def start_run(code: str, name: str, stopping_power_mode: int, force: bool) -> None:
	""" arrange the inputs for a LILAC or IRIS run in a machine-readable format and submit it to slurm
	    :param code: one of "LILAC" or "IRIS"
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	    :param stopping_power_mode: one of 0 for no stopping, 1 for Li-Petrasso-Zylstra, or 2 for Maynard-Deutsch
	    :param force: whether to run this even if another version of this run already exists TODO: let the user hit "y" so they don't have to redo the command
	"""
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
			run(["scancel", outputs_table.loc[name, "slurm ID"]])
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
	submission = run(["sbatch", script_path],
	                 check=True, capture_output=True, text=True)
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
	           np.stack([beam_radius, beam_intensity], axis=1),
	           delimiter=" ", fmt="%.6f")  # type: ignore

	# parse the materials
	shell_material = get_solid_material_from_name(inputs["shell material"])
	fill_material = get_gas_material_from_components(parse_gas_components(inputs["fill"]))

	pulse_end_time = pulse_time[pulse_power >= 1/2*np.max(pulse_power)][-1]

	# compose the input deck
	input_deck = build_lilac_input_deck(name, inputs, pulse_end_time, fill_material, shell_material)
	with open(f"{directory}/lilac_data_input.txt", "w") as file:
		file.write(input_deck)

	# compose the bash script
	bash_script = build_lilac_bash_script(name)
	with open(f"{directory}/run.sh", "w") as file:
		file.write(bash_script)

	return f"{directory}/run.sh"


def build_lilac_input_deck(
		name: str, inputs: Series, laser_off_time: float,
		fill_material: Material, shell_material: Material) -> str:
	""" construct the input deck corresponding to the given inputs and return it as a str """
	return fill_in_template(
		"lilac_input_deck.txt",
		parameters={
			# header
			"submitted": datetime.today().isoformat(" "),
			# rhydro namelist
			"sanitized name": sub(r"[/\\{}<> ~#%&*?]", "-", name),
			"absorption fraction": f"{inputs['absorption fraction']:.4f}",
			"nonthermal model": "none" if inputs["flux limiter"] > 0 else "vgon",
			"flux limiter": f"{inputs['flux limiter']:.3f}" if inputs["flux limiter"] > 0 else "0",
			"laser off time": f"{laser_off_time:.3f}",
			# first prof namelist (fill)
			"fill material code": f"{fill_material.material_code:d}",
			"fill protium percentage": f"{fill_material.protium_fraction*100:.2f}",
			"fill tritium percentage": f"{fill_material.tritium_fraction*100:.2f}",
			"fill EOS option": f"{fill_material.eos:d}",
			"fill opacity option": f"{fill_material.opacity:d}",
			"fill ionization option": f"{fill_material.ionization:d}",
			"fill pressure": f"{fill_material.pressure:.3f}",
			"fill radius": f"{inputs['outer diameter']/2 - inputs['shell thickness']:.2f}",
			# second prof namelist (shell)
			"shell material code": f"{shell_material.material_code:d}",
			"shell protium percentage": f"{shell_material.protium_fraction*100:.2f}",
			"shell tritium percentage": f"{shell_material.tritium_fraction*100:.2f}",
			"shell EOS option": f"{shell_material.eos:d}",
			"shell opacity option": f"{shell_material.opacity:d}",
			"shell ionization option": f"{shell_material.ionization:d}",
			"shell density": f"{shell_material.density:.3f}" if shell_material.density is not None else "None",
			"shell density multiplier": f"{inputs['density multiplier']:.4f}" if notnull(inputs["density multiplier"]) else "1",
			"shell thickness": f"{inputs['shell thickness']:.2f}",
			# third prof namelist (aluminum)
			"aluminum thickness": f"{inputs['aluminum thickness']:.2f}",
		},
		flags={
			"density specified": shell_material.density is not None,
			"aluminum": inputs["aluminum thickness"] > 0,
		}
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
				nuclide = nuclide_symbol(round(file["target/component_nuclear_charge"][i, j]),
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
			"number of times": f"{time.size:d}",
			"fill material code": f"{material_codes[0]:d}",
			"fill hydrogen fraction": f"{nuclide_fractions['H'][0]:.6f}",
			"fill deuterium fraction": f"{nuclide_fractions['D'][0]:.6f}",
			"fill tritium fraction": f"{nuclide_fractions['T'][0]:.6f}",
			"fill helium3 fraction": f"{nuclide_fractions['³He'][0]:.6f}",
			"fill carbon fraction": f"{nuclide_fractions['¹²C'][0]:.6f}",
			"fill mean atomic mass": f"{mean_atomic_masses[0]*1.66054e-27:.8g}",  # kg
			"fill max ionization": f"{mean_atomic_numbers[0]:.8g}",
			"shell material code": f"{material_codes[1]:d}",
			"shell hydrogen fraction": f"{nuclide_fractions['H'][1]:.6f}",
			"shell deuterium fraction": f"{nuclide_fractions['D'][1]:.6f}",
			"shell tritium fraction": f"{nuclide_fractions['T'][1]:.6f}",
			"shell helium3 fraction": f"{nuclide_fractions['³He'][1]:.6f}",
			"shell carbon fraction": f"{nuclide_fractions['¹²C'][1]:.6f}",
			"shell mean atomic mass": f"{mean_atomic_masses[1]*1.66054e-27:.8g}",  # kg
			"shell max ionization": f"{mean_atomic_numbers[1]}",
			"stopping power model": f"{stopping_power_mode:d}",
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
	args = parser.parse_args(sys.argv[2:])

	start_run(code, args.name, args.stopping_mode, args.force)
