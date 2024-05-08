import os.path
import sys
from argparse import ArgumentParser
from collections import defaultdict
from os import path
from typing import Sequence

import h5py
import numpy as np
from astropy.units import meter, centimeter, second, kilogram, joule, coulomb, farad, kiloelectronvolt
from matplotlib import pyplot as plt, colors
from numpy import stack, tile, cumsum, float64, nonzero, pi, average, exp, log, argmin, argmax, zeros
from pandas import Timestamp
from scipy import integrate

from data_io import write_row_to_outputs_table
from material import isotope_symbol, nuclide_symbol
from utilities import apparent_brightness, width, gradient, create_pdf

plt.rcParams["font.size"] = 11


def postprocess_lilac_run(name: str, lilac_status: str) -> None:
	""" read the raw results of a successful LILAC simulation and compile them into a human-readable
	    PDF.  also save some key numbers like yield and burn-average ρR to run_outputs.csv.
	"""
	# start by updating the run_outputs table
	write_row_to_outputs_table({
		"name": name,
		"status": lilac_status + " LILAC",
		"status changed": Timestamp.now(),
	})

	# check if there's any output
	directory = f"runs/{name}/lilac"
	if not os.path.isfile(f"{directory}/output.h5"):
		print("There was no LILAC output to postprocess.")
		return

	# if there is, get to postprocessing!
	print("loading LILAC output...")
	with h5py.File(f"{directory}/output.h5") as solution:
		print("processing LILAC output...")
		num_layers = solution["target/material_name"].size
		num_zones = solution["target/zone/atomic_mass"].size
		thermal_transport_model = solution["target/tt_model"][()].decode()

		# find the node positions
		node_position = solution["node/boundary_position"][:, :]  # (μm)
		zone_position = (node_position[0:-1, :] + node_position[1:, :])/2
		zone_height = node_position[1:, :] - node_position[0:-1, :]

		# get the compositions of the different layers
		mass_density = solution["zone/mass_density"][:, :]  # (g/cm^3)
		species_density = defaultdict(lambda: zeros(zone_position.shape))
		layers = []
		for i in range(num_layers):
			inner_index = solution["target/first_zone"][i] - 2  # subtract 2 because LILAC indexes from 2
			outer_index = solution["target/last_zone"][i] - 1
			thickness = node_position[outer_index, 0] - node_position[inner_index, 0]
			if solution["target/material_id"][i] == 13:
				layer_name = "Coating"
			elif solution["target/material_id"][i] == 102:
				layer_name = "Ice"
			elif mass_density[inner_index, 0] < .1:
				layer_name = "Fill"
			else:
				layer_name = "Shell"
			layers.append(Layer(
				layer_name,
				thickness,
				inner_index,
				outer_index,
				solution["target/material_name"][i].decode(),
				mass_density[inner_index, 0],
				solution["target/component_symbol"][i, :],
				solution["target/component_abundance"][i, :],
				solution["target/component_atomic_weight"][i, :],
				np.round(solution["target/component_nuclear_charge"][i, :]).astype(int),
			))

			# calculate the specific ion densities
			for j in range(len(layers[i].component_abundances)):
				species = nuclide_symbol(layers[i].atomic_numbers[j], round(layers[i].mass_numbers[j]))
				specific_density = layers[i].component_abundances[j]/(
						solution["target/atomic_weight"][i]*1.66054e-24)  # (g^-1)
				species_density[species][inner_index:outer_index, :] += \
					mass_density[inner_index:outer_index, :]*specific_density  # (cm^-3)

		interface_position = node_position[[layer.outer_index for layer in layers], :]

		time = solution["time/absolute_time"][:].astype(float64)  # (ns)
		laser_power = solution["time/laser_power"][:]/1e12  # (TW)
		# normalize the time so t=0 is the start of the laser pulse
		start_index = argmax(laser_power > laser_power.max()*1e-2)
		time = time - time[start_index]

		# pull out the different cumulative yield arrays
		fusion_reactions = ["DT-n", "D3He-p", "DD-n"]
		# note that LILAC calls these arrays "cumulative" but they're literally just not cumulative.
		# I convert them to cumulative because I think that's the most numericly precise way to do it.
		zone_cumulative_yield = {
			"DD-n": cumsum(solution["zone/cumulative_ddn_yield"][:, :], axis=1),  # neutrons/zone
			"DT-n": cumsum(solution["zone/cumulative_dtn_yield"][:, :], axis=1),  # neutrons/zone
			"D3He-p": cumsum(solution["zone/cumulative_dhe3_yield"][:, :], axis=1),  # protons/zone
		}

		# differentiate and integrate yield in time and space, respectively
		total_cumulative_yield, zone_yield_rate, total_yield_rate, total_yield = {}, {}, {}, {}
		for reaction in zone_cumulative_yield.keys():
			# change to double precision
			zone_cumulative_yield[reaction] = zone_cumulative_yield[reaction].astype(float64)
			# sum them in space
			total_cumulative_yield[reaction] = np.sum(zone_cumulative_yield[reaction], axis=0)  # neutrons
			# convert them to yield rates
			zone_yield_rate[reaction] = gradient(zone_cumulative_yield[reaction], time, axis=1)  # neutrons/zone/ns
			# convert them to yield rates
			total_yield_rate[reaction] = gradient(total_cumulative_yield[reaction], time, axis=0)  # neutrons/ns
			# and extract the final time-integrated value
			total_yield[reaction] = total_cumulative_yield[reaction][-1]
		main_fusion_reaction = max(fusion_reactions, key=lambda reaction: total_yield[reaction])

		# also calculate knock-on yields
		knockon_reactions = ["ko-d (peak)", "ko-t (peak)", "ko-p"]
		for reaction in knockon_reactions:
			particle = reaction[3]
			species_areal_density = np.sum(
				species_density[particle]*zone_height, axis=0)*1e-4  # particle/cm^2
			if particle == "p":
				knockon_cross_section = .6876e-24  # cm^2
			elif particle == "d":
				knockon_cross_section = .6440e-24*.157  # cm^2 (see C. K. Li and al., Phys. Plasmas 8 (2001) 4902)
			elif particle == "t":
				knockon_cross_section = .9403e-24*.135  # cm^2 (see C. K. Li and al., Phys. Plasmas 8 (2001) 4902)
			total_yield_rate[reaction] = total_yield_rate["DT-n"]*knockon_cross_section*species_areal_density*1.1  # particle/ns
			total_yield[reaction] = integrate.trapz(x=time, y=total_yield_rate[reaction])

		main_charged_particle = max(["D3He-p", "DD-n", "ko-d (peak)", "ko-p"], key=lambda particle: total_yield[particle])

		# calculate bang-time
		bang_index, bang_time, burn_width = {}, {}, {}
		for reaction in total_yield_rate.keys():
			if np.any(total_yield_rate[reaction] > 0):
				bang_index[reaction] = argmax(total_yield_rate[reaction])
				bang_time[reaction] = time[bang_index[reaction]]
				burn_width[reaction] = width(time, total_yield_rate[reaction])*1e3  # (ps)

		# set up the ion temperature in a way that allows easy pcolormeshing
		node_time = tile(time, (num_zones + 1, 1))
		ion_temperature = solution["zone/ion_temperature"][:, :]  # (keV)
		intertime_ion_temperature = (ion_temperature[:, 1:] + ion_temperature[:, 0:-1])/2

		electron_temperature = solution["zone/electron_temperature"][:, :].astype(float64)  # (keV)
		electron_density = solution["zone/electron_density"][:, :].astype(float64)  # (cm^-3)
		ionization = solution["zone/average_z"][:, :].astype(float64)

		# calculate the x-ray emission averaged temperature
		xray_energy_bands = []
		average_mass_density, average_electron_temperature, average_ionization = {}, {}, {}
		for cutoff in [1., 5., 20.]:  # (keV)
			key = f"{cutoff:.0f}+ keV x-ray"
			weights = solution["zone/zone_volume"][:, :]*apparent_brightness(
				np.maximum(0, electron_density), np.maximum(1e-3, electron_temperature), cutoff)
			if np.any(weights > 0):
				xray_energy_bands.append(key)
				bang_index[key] = argmax(np.sum(weights, axis=0))
				bang_time[key] = time[bang_index[key]]
				average_electron_temperature[key] = average(electron_temperature, weights=weights)
				average_mass_density[key] = average(mass_density, weights=weights)
				average_ionization[key] = average(ionization, weights=weights)

		# calculate coupling and degeneracy parameters
		ne = electron_density*centimeter**-3
		Te = electron_temperature*kiloelectronvolt
		e = 1.60e-19*coulomb  # (we're using units for just this part)
		ɛ0 = 8.85e-12*farad/meter
		ħ = 1.054e-34*joule*second
		me = 9.11e-31*kilogram
		coulomb_energy = e**2/(4*pi*ɛ0)*(4/3*pi*ne)**(1/3)
		fermi_energy = ħ**2/(2*me)*(3*pi**2*ne)**(2/3)
		coupling = coulomb_energy/(Te + fermi_energy)
		coupling = coupling.decompose().value
		degeneracy = Te/fermi_energy
		degeneracy = degeneracy.decompose().value

		# calculate the time-resolved areal density by layer
		areal_density = {}
		for i, interface_name in enumerate(["fill", "shell"]):
			start = solution["target/first_zone"][i] - 2
			end = solution["target/last_zone"][i] - 1
			areal_density[interface_name] = np.sum(
				(mass_density*zone_height)[start:end, :], axis=0)*1e-4*1e3  # (mg/cm^2)
		areal_density["total"] = areal_density["fill"] + areal_density["shell"]

	# calculate the stopping-averaged coupling parameter
	average_coupling, average_degeneracy = {}, {}
	time_weight = total_yield_rate[main_charged_particle]
	space_weight = mass_density/electron_temperature*zone_height
	space_weight = space_weight/np.sum(space_weight, axis=0)
	weights = time_weight*space_weight
	if np.any(weights > 0):
		average_electron_temperature["stopping"] = average(electron_temperature, weights=weights)
		average_mass_density["stopping"] = average(mass_density, weights=weights)
		average_coupling["stopping"] = exp(average(log(coupling), weights=weights))
		average_degeneracy["stopping"] = exp(average(log(degeneracy), weights=weights))

	average_ion_temperature = {}
	average_areal_density = {}
	for reaction in fusion_reactions:
		if np.any(zone_yield_rate[reaction] > 0):
			average_ion_temperature[reaction] = average(
				ion_temperature, weights=zone_yield_rate[reaction])
			average_areal_density[reaction] = average(
				areal_density["total"], weights=total_yield_rate[reaction])

	# in the event there is no fusion at all, make sure we have some values to report
	if main_fusion_reaction not in bang_time:
		# use time of minimum volume as a backup to bang-time
		bang_index[main_fusion_reaction] = argmin(interface_position[0, :])
		bang_time[main_fusion_reaction] = time[bang_index[main_fusion_reaction]]
		# and snapshot values at that time as a backup to averages
		average_ion_temperature[main_fusion_reaction] = ion_temperature[0, bang_index[main_fusion_reaction]]
		average_areal_density[main_fusion_reaction] = areal_density["total"][bang_index[main_fusion_reaction]]

	# calculate convergence ratio (defining radius with the gas-shell interface)
	radius = interface_position[0, :]
	convergence_ratio = radius[0]/np.min(radius)

	# save the pulse shape and burn history as a CSV
	table = stack(
		[time] + [total_yield_rate[reaction] for reaction in fusion_reactions],
		axis=1)
	np.savetxt(
		f"{directory}/time_history.csv", table, delimiter=",",  # type: ignore
		header="Time (ns), DD-n (ns^-1), DT-n (ns^-1), D3He-p (ns^-1)")

	# plot the spacio-temporally resolved ion temperature and shell trajectory
	print("generating default plots...")
	mm = 1/25.4
	fig, (ax_top_left, ax_bottom) = plt.subplots(
		2, 1, sharex="all", gridspec_kw=dict(hspace=0, height_ratios=[2, 3]), figsize=(140*mm, 120*mm), facecolor="none")
	ax_top_rite = ax_top_left.twinx()
	ax_bottom.set_facecolor("k")
	mesh = ax_bottom.pcolormesh(node_time, node_position, intertime_ion_temperature,
	                            cmap="inferno", norm=colors.LogNorm(vmin=.02, vmax=20))
	for i in range(num_layers):
		ax_bottom.plot(time, interface_position[i, :], 'w-', linewidth=1)
	ax_bottom.set_xlim(0, min(np.max(time), bang_time[main_fusion_reaction] + 1))
	ax_bottom.set_ylim(0, 1.2*interface_position[-1, 0])
	ax_bottom.set_xlabel(f"Time (ns)")
	ax_bottom.set_ylabel(f"Radius (μm)")
	fig.colorbar(mesh, ax=ax_bottom, orientation="horizontal", pad=.25,
	             label=f"Ion temperature (keV)")

	# plot the pulse shape and burn history
	peak = np.max([total_yield_rate[reaction] for reaction in total_yield_rate.keys()])
	curves, labels = [], []
	curves.append(
		ax_top_left.plot(time, laser_power, "k--", zorder=2)[0])
	labels.append("Laser")
	burn_histories_to_plot = total_yield_rate.keys()
	burn_histories_to_plot = filter(lambda key: total_yield[key] > 0, burn_histories_to_plot)
	burn_histories_to_plot = sorted(burn_histories_to_plot, key=lambda key: total_yield[key], reverse=True)
	for reaction in burn_histories_to_plot[:3]:
		if np.any(total_yield_rate[reaction] > 0):
			curves.append(
				ax_top_rite.plot(time, total_yield_rate[reaction], zorder=3)[0])
			labels.append(reaction)
	ax_top_left.set_xlim(0, min(np.max(time), bang_time[main_fusion_reaction] + 1))
	ax_top_left.set_ylim(0, None)
	ax_top_left.set_ylabel("Power (TW)")
	ax_top_left.locator_params(steps=[1, 2, 5, 10])
	ax_top_left.grid(axis="x")
	ax_top_rite.set_yscale("log")
	ax_top_rite.set_ylim(peak*1.5e-4, peak*1.5)
	ax_top_rite.set_ylabel("Yield (ns$^{-1}$)")
	ax_top_rite.grid(axis="y")
	ax_top_rite.legend(curves, labels, framealpha=1, fancybox=False)
	fig.tight_layout()
	fig.savefig(f"{directory}/time_plot.eps")
	fig.savefig(f"{directory}/time_plot.png", dpi=150)

	# plot the density and temperature profiles
	fig, (ax_top, ax_T) = plt.subplots(
		2, 1, sharex="all", gridspec_kw=dict(hspace=0), figsize=(140*mm, 100*mm), facecolor="none")
	i = bang_index[main_fusion_reaction]
	ax_T.plot(zone_position[:, i], ion_temperature[:, i], "C1--")
	ax_T.set_ylabel("$T_\\mathrm{i}$ (keV)")
	ax_T.set_ylim(0, None)
	ax_T.grid()
	ax_T.locator_params(steps=[1, 2, 5, 10])
	ax_T.set_xlabel("Radius (μm)")
	ax_ρ = ax_T.twinx()
	ax_ρ.plot(zone_position[:, i], mass_density[:, i], "C3-")
	ax_ρ.set_ylabel("$ρ$ (g/cm$^{3}$)")
	ax_ρ.set_ylim(0, None)
	ax_ρ.locator_params(steps=[1, 2, 5, 10])

	# plot the coupling and degeneracy parameters
	ax_top.plot(zone_position[:, i], coupling[:, i], "C0-.", label=r"$\Gamma_\mathrm{e}$")
	ax_top.plot(zone_position[:, i], degeneracy[:, i], "C2--", label=r"$\Theta_\mathrm{e}$")
	ax_top.legend(framealpha=1, fancybox=False)
	ax_top.set_xlabel("Radius (μm)")
	peak_mass_density = np.max(mass_density[:, i])
	outer_index = min(
		node_position.shape[0] - 1,
		nonzero(mass_density[:, i] > peak_mass_density/20)[0][-1] + 1)
	shell_radius = node_position[outer_index, i]
	ax_top.set_xlim(0, shell_radius)
	ax_top.set_yscale("log")
	ax_top.set_ylim(1e-2, 1e+2)
	ax_top.grid()
	fig.tight_layout()
	fig.savefig(f"{directory}/space_plot.eps")
	fig.savefig(f"{directory}/space_plot.png", dpi=150)

	# create the summary PDF!
	print("generating PDF...")
	pdf = create_pdf(f"{name} LILAC summary")
	# list the layers
	pdf.set_font("Noto", "B", 16)
	pdf.write(10, f"Capsule composition (outer diameter {2*interface_position[-1, 0]:.1f} μm)")
	pdf.ln()
	pdf.set_font("Noto", "", 16)
	for layer in reversed(layers):
		component_descriptions = []
		for i in range(len(layer.component_abundances)):
			if len(layer.component_names[i]) == 0:
				break
			if layer.component_abundances[i] > 0:
				component_descriptions.append(
					f"{layer.component_abundances[i]:.0%} {isotope_symbol(layer.atomic_numbers[i], layer.mass_numbers[i])}")
		layer_description = f"{layer.name}: {layer.thickness:.1f} μm {layer.material_name} ({' + '.join(component_descriptions)}) at {layer.density:.2g} g/cm³"
		pdf.set_x(15)
		pdf.write(10, layer_description)
		pdf.ln()
	# include the plots
	pdf.image(f"{directory}/time_plot.png", 10, 80, 140, 120)
	pdf.image(f"{directory}/space_plot.png", 150, 80, 140, 100)
	pdf.add_page()
	# print out some general numbers
	pdf.set_font("Noto", "B", 16)
	pdf.write(10, f"Overview quantities")
	pdf.ln()
	pdf.set_font("Noto", "", 16)
	pdf.set_x(15)
	pdf.write(10, f"Laser energy: {integrate.trapz(x=time, y=laser_power):.3f} kJ")
	pdf.ln()
	pdf.set_x(15)
	pdf.write(10, f"Convergence ratio: {convergence_ratio:.1f}")
	pdf.ln()
	pdf.set_x(15)
	pdf.write(10, f"Thermal transport: {thermal_transport_model:s}")
	pdf.ln()
	pdf.ln()
	# make pages for varius burn-averaged quantities
	averaged_quantities = [
		("Yield", total_yield, ".3g", ""),
		("Bang-time", bang_time, ".3f", "ns"),
		("Burn-width", burn_width, ".0f", "ps"),
		("T_ion", average_ion_temperature, ".2f", "keV"),
		("T_elec", average_electron_temperature, ".2f", "keV"),
		("ρ", average_mass_density, ".3g", "g/cm³"),
		("Z_eff", average_ionization, ".3f", ""),
		("Coupling", average_coupling, ".3g", ""),
		("Degeneracy", average_degeneracy, ".3g", ""),
		("ρR", average_areal_density, ".1f", "mg/cm²"),
	]
	for weighting in fusion_reactions + knockon_reactions + xray_energy_bands + ["stopping"]:
		if weighting not in total_yield or total_yield[weighting] > 0:
			pdf.set_font("Noto", "B", 16)
			pdf.write(10, f"{weighting} averaged quantities")
			pdf.ln()
			pdf.set_font("Noto", "", 16)
			for label, values, foremat, units in averaged_quantities:
				if weighting in values:
					pdf.set_x(15)
					pdf.write(
						10, f"{label}: {format(values[weighting], foremat)} {units}")
					pdf.ln()
			pdf.ln()
	# save it!
	pdf.output(f"{directory}/summary {path.split(name)[-1]}.pdf")

	# update our records
	print("updating run_outputs.csv...")
	write_row_to_outputs_table({
		"name": name,
		"yield": total_yield[main_fusion_reaction],
		"bang-time": bang_time[main_fusion_reaction],
		"convergence ratio": convergence_ratio,
		"areal density": average_areal_density[main_fusion_reaction],
		"ion temperature": average_ion_temperature[main_fusion_reaction],
	})
	print("done!")


class Layer:
	def __init__(self, name: str, thickness: float, inner_index: Sequence[float], outer_index: Sequence[float],
	             material_name: str, density: float,
	             component_names: Sequence[str], component_abundances: Sequence[float],
	             mass_numbers: Sequence[float], atomic_numbers: Sequence[int]):
		""" a layer of an unimploded capsule
			:param name: a human-readable identifier
		    :param thickness: the thickness (or radius in the gas fill's case) in μm
		    :param inner_index: the index of the node that marks its inner interface
		    :param outer_index: the index of the node that marks its outer interface
		    :param material_name: the name of the material
		    :param density: the initial density of the material in kg/L
		    :param component_names: the symbol of each nuclide in this material
		    :param component_abundances: the atomic fraction of each nuclide in this material
		    :param mass_numbers: the atomic weight of each nuclide in this material in Da
		    :param atomic_numbers: the atomic number of each nuclide in this material
		"""
		self.name = name
		self.thickness = thickness
		self.inner_index = inner_index
		self.outer_index = outer_index
		self.material_name = material_name
		self.density = density
		self.component_names = component_names
		self.component_abundances = component_abundances
		self.mass_numbers = mass_numbers
		self.atomic_numbers = atomic_numbers


if __name__ == "__main__":
	parser = ArgumentParser(
		prog="postprocess_lilac_run.sh",
		description = "read the raw results of a successful LILAC simulation and compile them into a human-readable PDF")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	parser.add_argument(
		"--status", type=str, default="completed",
		help="the end state of the run; one of 'completed', 'failed', or 'timeout'")
	args = parser.parse_args()

	postprocess_lilac_run(args.name, args.status)
	sys.exit(0)
