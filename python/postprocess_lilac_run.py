from argparse import ArgumentParser

import h5py
import numpy as np
from matplotlib import pyplot as plt, colors
from numpy import stack, tile, diff, cumsum, float64
from pandas import Timestamp

from python.data_io import load_outputs_table
from python.utilities import apparent_brightness, width, gradient


plt.rcParams["font.size"] = 12


def postprocess_lilac_run(name: str) -> None:
	""" read the raw results of a successful LILAC simulation and compile them into a human-readable
	    PDF.  also save some key numbers like yield and burn-average ρR to run_outputs.csv.
	"""
	directory = f"runs/{name}/lilac"
	with h5py.File(f'{directory}/output.h5') as solution:
		num_layers = solution["target/material_name"].size
		num_zones = solution["target/zone/atomic_mass"].size

		# start by getting the compositions of the different layers
		symbols = solution["target/component_symbol"][:, :]
		masses = solution["target/component_atomic_weight"][:, :]
		fractions = solution["target/component_abundance"][:, :]
		for i in range(num_layers):
			print(f"Layer {i}")
			j = 0
			while len(symbols[i, j]) > 0:
				print(f"  {symbols[i, j].decode()}-{masses[i, j]:.0f}: {fractions[i, j]:.4%}")
				j += 1
			print()

		time = solution["time/absolute_time"][:].astype(float64)  # (ns)
		laser_power = solution["time/laser_power"][:]/1e12  # (TW)
		# normalize the time so t=0 is the start of the laser pulse
		start_index = np.argmax(laser_power > laser_power.max()*1e-2)
		time = time - time[start_index]

		node_position = solution["node/boundary_position"][:, :]  # (μm)
		radius = node_position[-1, :]

		# pull out the different cumulative yield arrays
		reactions = {"DD-n", "DT-n", "D3He-p"}
		# note that LILAC calls these arrays "cumulative" but they're literally just not cumulative.
		zone_yield = {
			"DD-n": cumsum(solution["zone/cumulative_ddn_yield"][:, :], axis=1),  # neutrons/zone
			"DT-n": cumsum(solution["zone/cumulative_dtn_yield"][:, :], axis=1),  # neutrons/zone
			"D3He-p": cumsum(solution["zone/cumulative_dhe3_yield"][:, :], axis=1),  # protons/zone
		}
		# including total fusion yield
		zone_yield["all fusion"] = zone_yield["DD-n"] + zone_yield["DT-n"] + zone_yield["D3He-p"]

		# differentiate and integrate yield in time and space, respectively
		total_yield, zone_yield_rate, total_yield_rate = {}, {}, {}
		bang_time, burn_width = {}, {}
		for reaction in zone_yield.keys():
			# change to double precision
			zone_yield[reaction] = zone_yield[reaction].astype(float64)
			# sum them in space
			total_yield[reaction] = np.sum(zone_yield[reaction], axis=0)  # neutrons
			# convert them to yield rates
			zone_yield_rate[reaction] = gradient(zone_yield[reaction], time, axis=1)  # neutrons/zone/ns
			# convert them to yield rates
			total_yield_rate[reaction] = gradient(total_yield[reaction], time, axis=0)  # neutrons/ns
			# calculate bang-time
			if np.any(total_yield_rate[reaction] > 0):
				peak_emission = np.argmax(total_yield_rate[reaction])
				bang_time[reaction] = time[peak_emission]
				burn_width[reaction] = width(time, total_yield_rate[reaction])
		main_reaction = max(reactions, key=lambda reaction: total_yield[reaction][-1])

		# also calculate knock-on deuteron yield
		deuterium_areal_density = np.sum(
			solution["zone/deuteron_density"]*diff(node_position, axis=0), axis=0)*1e-4  # deuteron/cm^2
		knockon_cross_section = .100e-24  # cm^2
		total_yield_rate["ko-d"] = total_yield_rate["DT-n"]*deuterium_areal_density*knockon_cross_section  # deuteron/ns

		# set up the ion temperature in a way that allows easy pcolormeshing
		node_time = tile(time, (num_zones + 1, 1))
		ion_temperature = solution["zone/ion_temperature"][:, :]  # (keV)
		intertime_ion_temperature = (ion_temperature[:, 1:] + ion_temperature[:, 0:-1])/2

		plt.figure(facecolor="none", figsize=(8, 4))
		plt.gca().set_facecolor("k")
		plt.pcolormesh(node_time, node_position, intertime_ion_temperature,
		               cmap="inferno", norm=colors.LogNorm(vmin=.05, vmax=50))
		plt.xlim(0, min(np.max(time), bang_time["all fusion"] + 1))
		plt.ylim(0, 1.2*radius[0])
		plt.xlabel(f"Time (ns)")
		plt.ylabel(f"Position (μm)")
		plt.colorbar(label=f"Ion temperature (keV)")
		plt.title(name)
		plt.tight_layout()

		interface_positions = []
		for i in range(num_layers):
			interface_index = solution["target/last_zone"][i] - 1  # subtract 1 because LILAC indexes from 2
			interface_positions.append(node_position[interface_index, :])
			plt.plot(time, interface_positions[-1], 'w-', linewidth=1)

		plt.savefig(f'{directory}/rt_plot.png', dpi=300)

		filter_stacks = [
			[(50, "Al")],
			[(50, "Al"), (1, "SRIP"), (200, "Al")]
		]

		brightness = {}
		for filter_stack in filter_stacks:
			brightness[f"{filter_stack[0][0]}μm {filter_stack[0][1]} x-ray"] = apparent_brightness(
				solution["zone/average_z"][:, :],
				solution["zone/electron_density"][:, :],
				solution["zone/electron_temperature"][:, :],
				filter_stack)

		np.savetxt(f"{directory}/burn_history.csv",
		           stack([time] + [total_yield_rate[reaction] for reaction in reactions], axis=1),
		           delimiter=",",
		           header="Time (ns), DD-n (ns^-1), DT-n (ns^-1), D3He-p (ns^-1)")

		# create a laser and burn time history plot
		peak = np.max([total_yield_rate[reaction] for reaction in total_yield_rate.keys()])
		fig, ax_left = plt.subplots(facecolor="none", figsize=(8, 4))
		ax_left.plot(time, laser_power, "k--", label="Laser")
		ax_right = ax_left.twinx()
		for reaction in total_yield_rate.keys():
			if reaction != "all fusion":
				if np.any(total_yield_rate[reaction] > 0):
					ax_right.plot(time, total_yield_rate[reaction], label=reaction)
		ax_left.set_xlim(0, min(np.max(time), bang_time["all fusion"] + 1))
		ax_left.set_xlabel("Time (ns)")
		ax_left.set_ylim(0, None)
		ax_left.set_ylabel("Power (TW)")
		ax_right.set_yscale("log")
		ax_right.set_ylim(peak*1.5e-4, peak*1.5)
		ax_right.set_ylabel("Yield (ns^-1)")
		plt.tight_layout()
		plt.legend()
		plt.savefig(f"{directory}/burn.png", dpi=300)

		# calculate the time-resolved areal density by zone
		mass_density = solution["zone/mass_density"]
		areal_density = {}
		for i, name in enumerate(["fill", "shell"]):
			start = solution["target/first_zone"][i] - 2
			end = solution["target/last_zone"][i] - 1
			areal_density[name] = np.sum(
				(mass_density*diff(node_position, axis=0))[start:end, :], axis=0)*1e-4  # (g/cm^2)
		areal_density["total"] = areal_density["fill"] + areal_density["shell"]

		average_temperature = {}
		average_areal_density = {}
		for reaction in reactions:
			weights = zone_yield_rate[reaction]
			if np.all(weights == 0):
				continue
			temperature = solution["zone/ion_temperature"]
			average_temperature[reaction] = np.sum(temperature*weights)/np.sum(weights)
			weights = total_yield_rate[reaction]
			average_areal_density[reaction] = np.sum(areal_density["total"]*weights)/np.sum(weights)

	# calculate convergence ratio
	convergence_ratio = interface_positions[0][0]/np.min(interface_positions[0])

	# update our records
	outputs_table = load_outputs_table()
	outputs_table.loc[(name, "LILAC"), "status"] = "completed"
	outputs_table.loc[(name, "LILAC"), "status changed"] = Timestamp.now()
	outputs_table.loc[(name, "LILAC"), "yield"] = total_yield[main_reaction][-1]
	outputs_table.loc[(name, "LILAC"), "bang-time"] = bang_time[main_reaction]
	outputs_table.loc[(name, "LILAC"), "convergence ratio"] = convergence_ratio
	outputs_table.loc[(name, "LILAC"), "areal density"] = average_areal_density[main_reaction]
	outputs_table.loc[(name, "LILAC"), "ion temperature"] = average_temperature[main_reaction]
	outputs_table.to_csv("run_outputs.csv")


if __name__ == "__main__":
	parser = ArgumentParser(
		prog="postprocess_lilac_run.sh",
		description = "read the raw results of a successful LILAC simulation and compile them into a human-readable PDF")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	args = parser.parse_args()

	postprocess_lilac_run(args.name)
