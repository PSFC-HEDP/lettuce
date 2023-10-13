from __future__ import annotations

import json
from datetime import datetime
from re import fullmatch, sub, DOTALL, search
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from pandas import read_csv, DataFrame, Series

INPUT_DTYPES = {
	"name": str,
	"laser energy": float, "pulse shape": str, "beam profile": str,
	"outer diameter": float, "fill": str,
	"shell material": str, "shell thickness": float, "aluminum thickness": float,
	"absorption fraction": float, "flux limiter": float,
	"laser degradation": float, "density multiplier": float,
}
OUTPUT_DTYPES = {
	"name": str, "code": str,
	"status": str, "status changed": str, "slurm ID": str,
	"yield": float, "bang-time": float, "convergence ratio": float, "ρR": float,
	"ion temperature": float, "electron temperature": float,
}


LILAC_MATERIAL_CODES = {
	"H": 1, "He": 2, "Li": 3, "Be": 4, "C": 6, "diamond": 6, "HDC": 6,
	"Al": 13, "Si": 14, "Fe": 26, "Ge": 32, "Ta": 73, "Au": 79, "U": 92,
	"CH": 110, "CH2": 111, "strong CD": 113, "CD": 115, "CD2": 116,
	"SiO2": 150, "glass": 150, "polystyrene": 510,
}

LILAC_D3He_MIXTURES = {
	103: 1.00000, 289: 0.99000, 278: 0.92308, 292: 0.90000,
	280: 0.87500, 298: 0.80000, 291: 0.70000, 297: 0.66667,
	296: 0.60000, 295: 0.50000, 270: 0.45000, 283: 0.42857,
	288: 0.40000, 290: 0.33333, 279: 0.28571, 281: 0.28000,
	284: 0.20000, 282: 0.15000, 286: 0.11111, 287: 0.10000,
	285: 0.07692, 277: 0.02500, 299: 0.00990, 101: 0.00000,
}


def load_inputs_table() -> DataFrame:
	""" load the table containing the specifications for all of the runs """
	return read_csv(
		"run_inputs.csv", skipinitialspace=True, index_col="name", dtype=INPUT_DTYPES)


def load_outputs_table() -> DataFrame:
	""" load the table into which we will dump all of the run outputs """
	try:
		return read_csv(
			"run_outputs.csv", skipinitialspace=True, index_col=["name", "code"],
			dtype=OUTPUT_DTYPES, parse_dates=["status changed"])
	except IOError:
		table = DataFrame({key: Series(dtype=dtype) for key, dtype in OUTPUT_DTYPES.items()})
		table.set_index(["name", "code"], inplace=True)
		return table


def log_message(message: str) -> None:
	""" write a timestamped message to the runs log file, and also stdout """
	with open("runs.log", "a") as file:
		file.write(datetime.today().strftime('%m-%d %H:%M') + " | " + message + "\n")
	print(message)


def fill_in_template(template_filename: str, parameters: dict[str, str], flags: Optional[dict[str, bool]] = None) -> str:
	""" load a template from resources/templates and replace all of the angle-bracket-marked
	    parameter names with actual user-specified parameters
	    :param template_filename: the filename of the template to load, excluding resources/templates/
	    :param parameters: the set of strings that should be inserted into the <<>> spots
	    :param flags: a set of booleans to be used to evaluate <<if>> blocks
	    :return: a string containing the contents of the template with all the <<>>s evaluated
	    :raise KeyError: if there is a <<>> expression in the template and the corresponding value is not
	                     given in parameters or flags
	"""
	with open(f"resources/templates/{template_filename}") as template_file:
		content = template_file.read()

	# start with the flags
	if flags is not None:
		for key, active in flags.items():
			# if it's active, remove the <<if>> and <<endif>> lines
			if active:
				content = sub(f"<<if {key}>>\n", "", content)
				content = sub(f"<<endif {key}>>\n", "", content)
			# if it's inactive, remove the <<if>> and <<endif>> lines and everything between them
			else:
				content = sub(f"<<if {key}>>.*<<endif {key}>>\n", "", content, flags=DOTALL)
	# then do the parameter values
	for key, value in parameters.items():
		content = sub(f"<<{key}>>", value, content)

	# check to make sure we got it all
	remaining_blank = search("<<[a-z ]+>>", content)
	if remaining_blank:
		raise KeyError(f"you tried to fill out the template {template_filename} without specifying "
		               f"the value of {remaining_blank.group()}")

	return content


def load_pulse_shape(pulse_shape_name: str, total_energy: float) -> tuple[NDArray[float], NDArray[float]]:
	""" load a pulse shape from disk
	    :return: the time (ns) and total laser power (TW)
	"""
	filepath = f"resources/pulse_shapes/{pulse_shape_name}.json"
	try:
		with open(filepath, "r") as file:
			data = json.load(file)
	except IOError:
		raise IOError(f"the pulse shape file '{filepath}' is missing.  please "
		              f"download it from the OmegaOps pulse shape library.")
	if len(data) != 1:
		raise ValueError(f"the file '{filepath}' seems to contain multiple pulse "
		                 f"shapes.  please, when you download pulse shapes from OmegaOps, make sure you only have one "
		                 f"pusle shape activated each time you download a file.")
	if data[0]["Pulse"] != pulse_shape_name:
		raise ValueError(f"the file '{filepath}' seems to contain the information "
		                 f"for pulse shape {data[0]['Pulse']} instead of {pulse_shape_name}.  Please rename it "
		                 f"accordingly and download the true pulse shape file for {pulse_shape_name} from OmegaOps.")
	time = np.empty(len(data[0]["UV"]["data"]))
	power = np.empty(len(data[0]["UV"]["data"]))
	for i in range(len(data[0]["UV"]["data"])):
		time[i] = data[0]["UV"]["data"][i]["x"]/1e-9  # (convert to ns)
		power[i] = data[0]["UV"]["data"][i]["y"]
	raw_total = np.sum(power*np.gradient(time))
	power = power*total_energy/raw_total
	return time, power


def load_beam_profile(beam_profile_name: str) -> tuple[NDArray[float], NDArray[float]]:
	""" load a beam profile from disk
	    :return: the radial coordinate (μm) and laser beam intensity (normalized)
	"""
	filepath = f"resources/beam_profiles/{beam_profile_name.replace(' ', '_')}.txt"
	try:
		data = np.loadtxt(filepath)
	except IOError:
		raise IOError(f"the beam profile file '{filepath}' is missing.  please ask Varchas for appropriate data.")
	if data.shape[1] != 2 or data[-1, 0] != -1 or data[-1, 1] != 0:
		raise ValueError(f"the file '{filepath}' seems wrong.  it should be a two-column whitespace-separated value "
		                 f"file with '-1 0' for the last row.")
	return data[:, 0], data[:, 1]


def find_best_D3He_material_code(fHe: float) -> tuple[int, float]:
	""" find the LILAC material code that represents the mixture of D and 3He that most accurately
	    represents the specified fill ratio.  specifically, find the LILAC material with the nearest
	    atomic 3He fraction that is less than or equal to the given one.
	    :return: the selected material code and the atomic 3He fraction that corresponds to it
	"""
	if 0 <= fHe <= 1:
		for matcode in LILAC_D3He_MIXTURES:
			if LILAC_D3He_MIXTURES[matcode] <= fHe + 1e-3:
				return matcode, LILAC_D3He_MIXTURES[matcode]
	raise ValueError(f"something's wrong with the 3He fill fraction {fHe:.0%}")


def get_shell_material_from_name(name: str) -> Material:
	""" create the Material object that describes the solid named by name.  see the [LILAC material
	    table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for
	    more information.
	"""
	if name == "D" or name == "DD" or name == "D2":
		return Material(102)
	elif name == "H" or name == "HH" or name == "H2":
		return Material(102, protium_fraction=1.0)
	elif name == "T" or name == "TT" or name == "T2":
		return Material(102, tritium_fraction=1.0)
	elif name == "DT":
		return Material(102, tritium_fraction=0.5)
	elif name in LILAC_MATERIAL_CODES:
		return Material(LILAC_MATERIAL_CODES[name])
	else:
		try:
			return Material(int(name))
		except ValueError:
			raise IndexError(f"I don't recognize the material '{name}'")


def get_gas_material_from_components(partial_pressures: dict[str, float]) -> Material:
	""" create the Material object that describes the solid named by name.  see the [LILAC material
	    table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for
	    more information.
	    :param partial_pressures: a dictionary where each key is the name of a nuclide and each
	                              value is its partial pressure in atm
	"""
	total_pressure = sum(partial_pressures.values())
	nuclides = partial_pressures.keys()

	# convert the molecular partial pressures to atomic fractions
	atomic_fractions = {}
	for nuclide in nuclides:
		if nuclide in ["H", "D", "T", "N", "O", "F", "Cl"]:
			atomic_fractions[nuclide] = 2*partial_pressures[nuclide]
		else:
			atomic_fractions[nuclide] = partial_pressures[nuclide]
	normalization = sum(atomic_fractions.values())
	for nuclide in nuclides:
		atomic_fractions[nuclide] /= normalization

	# the way gasses are specified is an uscritic mess, so let me explain...
	# if this is just a mixture of hydrogen isotopes
	if all(nuclide in ["H", "D", "T"] for nuclide in nuclides):
		# we can use bilt-in gas #101 and control ratios using ptrit
		return Material(101,
		                protium_fraction=atomic_fractions.get("H", 0.),
		                tritium_fraction=atomic_fractions.get("T", 0.),
		                pressure=total_pressure)
	# if this is a mixture of hydrogen and helium-3
	elif all(nuclide in ["H", "D", "T", "3He"] for nuclide in nuclides):
		# we can use a bilt-in gas iff there's a matcode that exactly matches the desired ratio
		best_material_code, found_atomic_fraction = find_best_D3He_material_code(atomic_fractions["3He"])
		if abs(found_atomic_fraction - atomic_fractions["3He"]) < 1e-2:
			return Material(best_material_code,
			                protium_fraction=atomic_fractions.get("H", 0.),
			                tritium_fraction=atomic_fractions.get("T", 0.),
			                pressure=total_pressure)
		else:
			# otherwise, we will need to specify this material manually using `mater`
			print("warning: I don't think this method of specifying the fill fraction, like, actually works...")
			return Material(-best_material_code,
			                protium_fraction=atomic_fractions.get("H", 0.),
			                tritium_fraction=atomic_fractions.get("T", 0.),
			                pressure=total_pressure)

	else:
		raise NotImplementedError("I don't autodetect materials by species unless it's a simple gas mix.")


def parse_gas_components(descriptor: str) -> dict[str, float]:
	""" read a string that contains information about the components and pressure of a gas mixture
	    :param descriptor: a string containing a list of nuclides and partial pressures.  it should look something
	                       like this: "12atm 3He + 6atm D". Each component is given as a molecular pressure at 293K
	                       followed by the name of the element. Note that these are molecular pressures. "D" and "D2"
	                       are interchangeable.
	    :raise ValueError: if the given descriptor cannot be parsed for whatever reason
	"""
	parts = descriptor.split("+")
	components = {}
	for part in parts:
		reading = fullmatch(r"\s*([0-9.]+)(\s?atm)?\s*([0-9]*[A-Z][a-z]?)2?\s*", part)
		if reading is None:
			raise ValueError(f"cannot parse '{part}'")
		pressure = float(reading.group(1))
		nuclide = reading.group(3)
		if nuclide in components:
			raise ValueError(f"the nuclide '{nuclide}' seems to appear twice in '{descriptor}'.")
		components[nuclide] = pressure
	return components


class Material:
	def __init__(
			self, material_code: int,
			protium_fraction: float = 0.0, tritium_fraction: float = 0.0,
			density: Optional[float] = None, pressure: Optional[float] = None):
		""" define a material by its composition and density
		    :param material_code: the LILAC material code (see LILAC user guide)
		    :param protium_fraction: the fraction of atoms that are protium.  this must be less
		                             than or equal to the atomic deuterium fraction.
		    :param tritium_fraction: the fraction of atoms that are tritium.  this must be less
		                             than or equal to the atomic deuterium fraction.
		    :param density: the density in g/mL
		    :param pressure: the pressure at room-temperature in Pa
		"""
		self.material_code = material_code
		self.protium_fraction = protium_fraction
		self.tritium_fraction = tritium_fraction
		self.density = density
		self.pressure = pressure

	def __str__(self):
		return f"Material({self.material_code}, protium_fraction={self.protium_fraction}, " \
		       f"tritium_fraction={self.tritium_fraction}, density={self.density}, " \
		       f"pressure={self.pressure})"

	def __eq__(self, other: Material):
		return self.material_code == other.material_code and \
		       self.protium_fraction == other.protium_fraction and \
		       self.tritium_fraction == other.tritium_fraction and \
		       self.density == other.density and \
		       self.pressure == other.pressure
