import os.path
import sys
from argparse import ArgumentParser
from os import path

import h5py
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from numpy import diff, ndarray, repeat
from pandas import Timestamp

from data_io import write_row_to_outputs_table
from utilities import create_pdf

plt.rcParams["font.size"] = 11


def postprocess_iris_run(name: str, iris_status: str) -> None:
	""" read the raw results of a successful IRIS simulation and compile them into a human-readable
	    PDF with plots of the image and spectrum.
	"""
	# start by updating the run_outputs table
	write_row_to_outputs_table({
		"name": name,
		"status": iris_status + " IRIS",
		"status changed": Timestamp.now(),
	})

	# check if there's any output
	directory = f"runs/{name}/iris-model1"
	if not os.path.isfile(f"{directory}/output.h5"):
		print("There was no IRIS output to postprocess.")
		return

	# if there is, get to postprocessing!
	print("loading IRIS output...")
	with h5py.File(f"{directory}/output.h5") as solution:
		image = solution["images/image"][:, 0:6, :, 0, :, :].sum(axis=(0, 1, 2))  # take a single imager, and sum all the neutron images over time, and energy
		spacial_range = solution["images/size"][0]/1e-6  # and take the size of the image in μm
		DT_n_spectrum = solution["spectra/dNdE"][:, 0, :, 0:6, :, 0].sum(axis=(0, 1, 2))  # take a single spectrometer, and sum all neutron types over all time
		DD_n_spectrum = solution["spectra/dNdE"][:, 1, :, 0:6, :, 0].sum(axis=(0, 1, 2))  # take a single spectrometer, and sum all neutron types over all time
		ko_d_spectrum = solution["spectra/dNdE"][:, :, :, 6, :, 0].sum(axis=(0, 1, 2))  # take a single spectrometer, and sum the deuteron spectrum over all time
		ko_t_spectrum = solution["spectra/dNdE"][:, :, :, 7, :, 0].sum(axis=(0, 1, 2))  # take a single spectrometer, and sum all neutron types over all time
		energy_bins = solution["spectra/energy"][:]  # and take the energy bins in MeV

	# plot the neutron image and possibly also a charged particle image
	print("generating default plots...")
	mm = 1/25.4
	fig = plt.figure(figsize=(140*mm, 130*mm), facecolor="none")
	ax = fig.add_subplot()
	ax.imshow(image.T, origin="lower",
	          extent=(-spacial_range, spacial_range, -spacial_range, spacial_range),
	          cmap="magma")
	ax.set_xlabel(f"x (μm)")
	ax.set_ylabel(f"y (μm)")
	fig.tight_layout()
	fig.savefig(f"{directory}/image_plot.eps")
	fig.savefig(f"{directory}/image_plot.png", dpi=150)

	# plot the neutron, deuteron, and triton spectra
	fig = plt.figure(figsize=(140*mm, 130*mm), facecolor="none")
	ax = fig.add_subplot()
	ax.grid()
	plot_histogram(ax, energy_bins/1.60e-13, DD_n_spectrum + DT_n_spectrum, alpha=1.0, label="Neutrons")
	if np.sum(ko_d_spectrum) > 0:
		plot_histogram(ax, energy_bins/1.60e-13, ko_d_spectrum, alpha=.7, label="Deuterons")
	if np.sum(ko_t_spectrum) > 0:
		plot_histogram(ax, energy_bins/1.60e-13, ko_t_spectrum, alpha=.7, label="Tritons")
	ax.set_xlabel("Energy (MeV)")
	ax.set_xlim(0, 16)
	ax.set_yscale("log")
	peak_magnitude = max(np.max(DD_n_spectrum + DT_n_spectrum), np.max(ko_d_spectrum), np.max(ko_t_spectrum))
	ax.set_ylim(peak_magnitude*2e-6, peak_magnitude*2e+0)
	fig.tight_layout()
	ax.legend(loc="upper left")
	fig.savefig(f"{directory}/spectrum_plot.eps")
	fig.savefig(f"{directory}/spectrum_plot.png", dpi=150)

	# create the summary PDF!
	print("generating PDF...")
	pdf = create_pdf(f"{name} IRIS summary")
	# include the plots
	pdf.image(f"{directory}/image_plot.png", 10, 20, 140, 130)
	pdf.image(f"{directory}/spectrum_plot.png", 150, 20, 140, 130)
	# print out some yields
	pdf.set_y(150)
	pdf.set_font("Noto", "B", 16)
	pdf.write(10, f"Yields")
	pdf.ln()
	pdf.set_font("Noto", "", 16)
	pdf.set_x(15)
	pdf.write(10, f"DT-neutron: {np.sum(DT_n_spectrum*diff(energy_bins)):.4g} ")
	pdf.ln()
	pdf.set_x(15)
	pdf.write(10, f"DD-neutron: {np.sum(DT_n_spectrum*diff(energy_bins)):.4g}")
	pdf.ln()
	pdf.set_x(15)
	pdf.write(10, f"ko-deuteron: {np.sum(ko_d_spectrum*diff(energy_bins)):.4g}")
	pdf.ln()
	pdf.set_x(15)
	pdf.write(10, f"ko-triton: {np.sum(ko_t_spectrum*diff(energy_bins)):.4g}")
	pdf.ln()
	# save it!
	pdf.output(f"{directory}/summary {path.split(name)[-1]}.pdf")

	print("done!")


def plot_histogram(ax: Axes, bin_edges: ndarray, densities: ndarray, **kwargs):
	ax.fill_between(repeat(bin_edges, 2)[1:-1], 0, repeat(densities, 2),
	                edgecolor="none", **kwargs)


if __name__ == "__main__":
	parser = ArgumentParser(
		prog="postprocess_iris_run.sh",
		description = "read the raw results of a successful IRIS simulation and compile them into a human-readable PDF")
	parser.add_argument(
		"name", type=str,
		help="the name of the run, as specified in run_inputs.csv")
	parser.add_argument(
		"--status", type=str, default="completed",
		help="the end state of the run; one of 'completed', 'failed', or 'timeout'")
	args = parser.parse_args()

	postprocess_iris_run(args.name, args.status)
	sys.exit(0)
