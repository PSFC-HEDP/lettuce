from __future__ import annotations

import json
import re
from typing import Optional

import numpy as np
import pandas as pd
import periodictable as periodictable


D_3He_MIXTURES = {
	103: 1.00000, 289: 0.99000, 278: 0.92308, 292: 0.90000,
	280: 0.87500, 298: 0.80000, 291: 0.70000, 297: 0.66667,
	296: 0.60000, 295: 0.50000, 270: 0.45000, 283: 0.42857,
	288: 0.40000, 290: 0.33333, 279: 0.28571, 281: 0.28000,
	284: 0.20000, 282: 0.15000, 286: 0.11111, 287: 0.10000,
	285: 0.07692, 277: 0.02500, 299: 0.00990, 101: 0.00000,
}


def load_pulse_shape(pulse_shape_name, total_energy):
	""" load a pulse shape from disk """
	try:
		with open(f"pulse_shapes/{pulse_shape_name}.json", "r") as file:
			data = json.load(file)
	except IOError:
		raise ValueError(F"")
	if len(data) != 1:
		raise ValueError(f"the file 'pulse_shapes/{pulse_shape_name}.json' seems to contain multiple pulse shapes.  please, when you download pulse shapes from OmegaOps, make sure you only have one pule shape activated each time you download a file.")
	if data[0]["Pulse"] != pulse_shape_name:
		raise ValueError(f"the file 'pulse_shapes/{pulse_shape_name}.json' seems to contain the information for pulse shape {data[0]['Pulse']} instead of {pulse_shape_name}.  Please rename it accordingly and download the true pulse shape file for {pulse_shape_name} from OmegaOps.")
	time = np.empty(len(data[0]["UV"]["data"]))
	power = np.empty(len(data[0]["UV"]["data"]))
	for i in range(len(data[0]["UV"]["data"])):
		time[i] = data[0]["UV"]["data"][i]["x"]/1e-9  # (convert to ns)
		power[i] = data[0]["UV"]["data"][i]["y"]
	raw_total = np.sum(power*np.gradient(time))
	power = power*total_energy/raw_total
	return time, power


def find_best_d3he_material_code(fHe):
	for matcode in D_3He_MIXTURES:
		if D_3He_MIXTURES[matcode] <= fHe + 1e-3:
			return matcode
	return ValueError(f"something's wrong with the 3He fill fraction {fHe:.0%}")


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
	elif name == "CH":
		return Material(110)
	elif name == "strong CD":
		return Material(113)
	elif name == "CD":
		return Material(115)
	elif name == "polystyrene":
		return Material(510)
	elif name == "SiO2" or name == "glass":
		return Material(150)
	else:
		try:
			# Symbol is given as carbon, hydrogen, etc.
			table_entry = periodictable.elements.symbol(name)
		except ValueError:
			raise IndexError(f"Couldnt find entry for '{name}'")
		else:
			return Material(table_entry.number)


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
		return Material(
			101, protium_fraction=atomic_fractions.get("H", 0),
			tritium_fraction=atomic_fractions.get("T", 0),
			pressure=total_pressure)
	# if this is a mixture of hydrogen and helium-3
	elif all(nuclide in ["H", "D", "T", "He3"] for nuclide in nuclides):
		# we can use a bilt-in gas iff there's a matcode that exactly matches the desired ratio
		best_material_code, found_atomic_fraction = find_best_d3he_material_code(atomic_fractions["He3"])
		if abs(found_atomic_fraction - atomic_fractions["3He"]) < 1e-2:
			return Material(best_material_code,
			                protium_fraction=atomic_fractions.get("H", 0),
			                tritium_fraction=atomic_fractions.get("T", 0),
			                pressure=total_pressure)
		else:
			# otherwise, we will need to specify this material manually using `mater`
			print("warning: I don't think this method of specifying the fill fraction, like, actually works...")
			return Material(-best_material_code,
			                protium_fraction=atomic_fractions.get("H", 0),
			                tritium_fraction=atomic_fractions.get("T", 0),
			                pressure=total_pressure)

	else:
		raise NotImplementedError("I don't autodetect materials by species unless it's a simple gas mix.  please specify the name or code. for custom gas mixen, use code -999.")


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
		reading = re.fullmatch(r"\s*([0-9.]+)(\s?atm)?\s*([0-9]*[A-Z][a-z]?)2?\s*", part)
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
