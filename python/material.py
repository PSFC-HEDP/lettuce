from __future__ import annotations

from copy import copy
from collections import defaultdict
from re import fullmatch
from typing import Optional, Union

from numpy import array, argmin

from utilities import to_superscript, from_superscript, drop_zeros


class Material:
	def __init__(
			self, material_code: int, *,
			eos: int, opacity: Union[int, str], ionization: int,
			density: Optional[float] = None, pressure: Optional[float] = None,
			components: Optional[dict[str, float]] = None,
			opacity_table: Optional[str] = None):
		""" define a material by its composition and density
		    :param material_code: the LILAC material code (see LILAC user guide).  it may be negative instead of
		                          positive to indicate that LILAC doesn't have an exact match, but the absolute value
		                          of the material code should be the closest material that LILAC does support.
		    :param eos: the equation of state table option number:
		                1 - ideal gas
		                2 - perfect degenerate gas
		                4 - SESAME
		                5 - LLNL table
		                6 - analytic thomas-fermi
		                7 - QEOS
		                8 - first principle eos
		                9 - 2011 LLNL table for D2
		    :param opacity: either the filename of the opacity table to use or the opacity table option number:
		                    0 - none
		                    1 - astrophysical table, planckian mean ('aot_pl_*.prp')
		                    2 - none but changed to 20 internally
		                    7 - astrophysical table, planckian mean, generated at runtime with aplmix
		                    8 - first principle table
		                    9 - LLNL table (only 101, 102, 110, 115)
		                    10 - astrophysical table, rosseland mean ('aot_ro_*.prp')
		                    20 - NLTE average-ion
		                    21 - tabular non-LTE ('avi_pl_*.prp')
		                    22 - tabular LTE ('ave_pl_*.prp')
		                    23 - tabular non-LTE from Reuben ('psi_pl_*.prp')
		                    27 - Prism LTE ('plt_*.prp')
		                    28 - Prism collisional radiative equilibrium non-LTE ('pcr_*.prp')
		                    29 - Native reader for Prism non-LTE ('prp_*.prp')
		                    40 - constant opacity
		    :param ionization: the ionization table option number:
		                       1 - astrophysical table
		                       2 - average-ion model
		                       3 - thomas-fermi table
		                       4 - tabular non-LTE
		                       5 - tabular LTE
		                       6 - Orchid analytic formulation
		                       8 - first-principle opacity table
		                       10 - fully ionized
		    :param density: the density in g/mL
		    :param pressure: the pressure at room-temperature in Pa
		    :param components: a dictionary specifying what ion species are present in the material and their atomic
		                       abundance fractions.  each key should be the atomic symbol of the ion's element; the mass
		                       will be assumed to be the natural average mass number.  a specific mass number may be
		                       prepended.  for custom materials the values should add to 1.  for materials LILAC already
		                       knows about, this argument may be omitted or passed with only the H and T abundances.
		                       yes, H and ¹H are technicly different, but the actual mass of a proton is so close to the
		                       atomic mass of hydrogen it doesn't matter.
		"""
		self.material_code = material_code
		self.density = density
		self.pressure = pressure
		self.eos = eos
		self.opacity = opacity
		self.ionization = ionization
		self.components = components if components is not None else {}
		# determine the default opacity table if it's not given
		if opacity_table is not None:
			self.opacity_table = opacity_table
		elif opacity == 1:
			self.opacity_table = f"aot_pl_48_50x50Dt_{abs(self.material_code)}.txt"
		elif opacity == 8:
			self.opacity_table = f"FPOT_pl_48_50x50Dt_{abs(self.material_code)}.txt"
		elif opacity == 10:
			self.opacity_table = f"aot_ro_48_50x50Dt_{abs(self.material_code)}.txt"
		elif opacity == 21:
			self.opacity_table = f"avi_pl_48_50x50Dt_{abs(self.material_code)}.txt"
		elif opacity == 22:
			self.opacity_table = f"ave_pl_48_50x50Dt_{abs(self.material_code)}.txt"
		else:  # if it's for an opacity option I (Justin) don't understand
			self.opacity_table = None  # just hope that it doesn't come up


	def __eq__(self, other: Material) -> bool:
		return self.material_code == other.material_code and \
		       self.components == other.components and \
		       self.density == other.density and \
		       self.pressure == other.pressure and \
		       self.eos == other.eos and \
		       self.opacity == other.opacity and \
		       self.ionization == other.ionization

	def __str__(self) -> str:
		return repr(self)

	def __repr__(self) -> str:
		return f"Material({self.material_code!r}, eos={self.eos!r}, opacity={self.opacity!r}, " \
		       f"ionization={self.ionization!r}, components={self.components!r}, " \
		       f"density={self.density!r}, pressure={self.pressure!r})"


# list of atomic symbols
ATOMIC_SYMBOLS = [
	"n", "H", "He",
	"Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
	"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
	"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
	"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
	"Ac", "Th", "Pa", "U", "Np", "Pu"]
ATOMIC_MASSES = [
	1.007, 1.008, 4.003,
	6.940, 9.012, 10.810, 12.011, 14.007, 15.999, 18.998, 20.180,
	22.990, 24.305, 26.982, 28.085, 30.974, 32.060, 35.450, 39.948,
	39.098, 40.078, 44.956, 47.867, 50.941, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.380, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798,
	85.468, 87.620, 88.906, 91.224, 92.906, 95.950, 97.907, 101.070, 102.906, 106.420, 107.868, 112.414, 114.818, 118.710, 121.760, 127.600, 126.904, 131.293,
	132.905, 137.327, 138.905, 140.116, 140.908, 144.242, 144.913, 150.360, 151.964, 157.250, 158.925, 162.500, 164.930, 167.259, 168.934, 173.045, 174.967, 178.490, 180.948, 183.840, 186.207, 190.230, 192.217, 195.084, 196.967, 200.592, 204.380, 207.200, 208.980, 209.000, 210.000, 222.000,
	223.000, 226.000, 227.000, 232.038, 231.036, 238.029,
]


# how to specify various common solids in LILAC
# see the [LILAC material table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for more information.
LILAC_SOLID_MATERIALS = {
	"H":           Material(102, eos=8, ionization=1, opacity=8, components={"¹H": 1.00},
	                        opacity_table="HDT_plt_50x50Dt_101.prp"),
	"DT":          Material(102, eos=8, ionization=1, opacity=8, components={"T": 0.50},
	                        opacity_table="HDT_plt_50x50Dt_101.prp"),
	"Li":          Material(3, eos=6, ionization=2, opacity=1),
	"Be":          Material(4, eos=6, ionization=2, opacity=1),
	"diamond":     Material(6, eos=6, ionization=2, opacity=1, density=3.50),
	"HDC":         Material(6, eos=6, ionization=2, opacity=1, density=3.50),
	"Al":          Material(13, eos=6, ionization=2, opacity=1,
	                        opacity_table="Al_pcr_50x50Dt_013.prp"),
	"Si":          Material(14, eos=6, ionization=2, opacity=1),
	"Fe":          Material(26, eos=6, ionization=2, opacity=1),
	"Ge":          Material(32, eos=6, ionization=2, opacity=1),
	"Ta":          Material(73, eos=6, ionization=2, opacity=1),
	"Au":          Material(79, eos=6, ionization=4, opacity=21),
	"U":           Material(92, eos=6, ionization=2, opacity=1),
	"CH":          Material(110, eos=8, ionization=1, opacity=8, density=1.03,
	                        opacity_table="CH_pcr_50x50Dt_110.prp"),
	"strong CD":   Material(113, eos=8, ionization=1, opacity=8, density=1.09,
	                        opacity_table="CD_pcr_50x50Dt_113.prp"),
	"CD":          Material(115, eos=8, ionization=1, opacity=1, density=1.10),
	"CHTi":        Material(125, eos=6, ionization=1, opacity=1),
	"CHGe":        Material(148, eos=8, ionization=2, opacity=21, density=1.11,
	                        opacity_table="CHGe(0p94%)_pcr_50x50DT_a.prp"),
	"SiO2":        Material(150, eos=6, ionization=1, opacity=1, density=2.20,
	                        opacity_table="SiO2_pcr_50x50Dt_150.prp"),
	"glass":       Material(150, eos=6, ionization=1, opacity=1, density=2.20),
	"CHCu":        Material(232, eos=8, ionization=1, opacity=1, density=1.23),
	"CHSi":        Material(356, eos=4, ionization=4, opacity=21, density=1.24,
	                        opacity_table="CHSi(7p4%)_pcr_50x50Dt_354.prp"),
	"polystyrene": Material(510, eos=8, ionization=1, opacity=8, density=1.05,
	                        opacity_table="CH_pcr_50x50Dt_510.prp"),
}

# how to specify various common gases in LILAC
LILAC_GAS_MATERIALS = {
	"D": Material(101, eos=8, ionization=1, opacity=8,
	              components={"D": 1}),
	"DT": Material(101, eos=8, ionization=1, opacity=8,
	               components={"D": .50, "T": .50},
	               opacity_table="HDT_plt_50x50Dt_101.prp"),
	"He": Material(2, eos=6, ionization=2, opacity=1),
	"N": Material(7, eos=6, ionization=2, opacity=1),
	"Ne": Material(10, eos=6, ionization=2, opacity=1),
	"Ar": Material(18, eos=8, ionization=1, opacity=20),
	"D³He": Material(295, eos=6, ionization=1, opacity=1,
	                 components={"D": .50, "³He": .50}),
	"D3He": Material(295, eos=6, ionization=1, opacity=1,
	                 components={"D": .50, "³He": .50}),
	"air": Material(300, eos=4, ionization=1, opacity=1,
	                components={"N": .78, "O": .21, "Ar": .01},
	                opacity_table="Air_pcr_50x50Dt_300.prp"),
}

# material codes that correspond to mixtures of D³He (values are atomic ³He fraction)
LILAC_D3He_MIXTURES = {
	103: 1.00000, 289: 0.99000, 278: 0.92308, 292: 0.90000,
	280: 0.87500, 298: 0.80000, 291: 0.70000, 297: 0.66667,
	296: 0.60000, 295: 0.50000, 283: 0.42857, 288: 0.40000,
	290: 0.33333, 279: 0.28571, 281: 0.28000, 284: 0.20000,
	282: 0.15000, 286: 0.11111, 287: 0.10000, 285: 0.07692,
	277: 0.02500, 299: 0.00990, 101: 0.00000,
}


def nuclide_symbol(atomic_number: int, mass_number: int) -> str:
	""" succinctly describe this particular nuclide (p, d, ³He, and so on) """
	if atomic_number == -1:
		return "e"
	elif atomic_number == 1 and mass_number == 1:
		return "p"
	elif atomic_number == 1 and mass_number == 2:
		return "d"
	elif atomic_number == 1 and mass_number == 3:
		return "t"
	else:
		return to_superscript(str(mass_number)) + ATOMIC_SYMBOLS[atomic_number]


def isotope_symbol(atomic_number: int, mass_number: float) -> str:
	""" succinctly describe this particular isotope (¹H, D, ³He, and so on) """
	if atomic_number == 1 and mass_number == 2:
		return "D"
	elif atomic_number == 1 and mass_number == 3:
		return "T"
	elif mass_number == ATOMIC_MASSES[atomic_number]:
		return ATOMIC_SYMBOLS[atomic_number]
	elif mass_number == int(mass_number):
		return to_superscript(str(int(mass_number))) + ATOMIC_SYMBOLS[atomic_number]
	else:
		return str(mass_number) + ATOMIC_SYMBOLS[atomic_number]


def parse_isotope_symbol(symbol: str) -> tuple[int, float]:
	""" take a succinct description of an isotope and deduce its atomic number and mass number
	    :return: the atomic number followed by the atomic weight in Da
	"""
	if symbol == "D":
		return 1, 2.014
	elif symbol == "T":
		return 1, 3.016
	else:
		reading = fullmatch(r"([0-9.]*)([A-Z][a-z]?)", from_superscript(symbol))
		if reading is None:
			raise ValueError(f"'{symbol}' does not appear to be a nuclide symbol.")
		mass_number, element_symbol = reading.groups()
		atomic_number = ATOMIC_SYMBOLS.index(element_symbol)
		if len(mass_number) > 0:
			mass_number = float(mass_number)  # take the mass number at face value (ignore the binding energy)
		else:
			mass_number = ATOMIC_MASSES[atomic_number]  # or if none is given, use the natural atomic mass
		return atomic_number, mass_number


def expand_compound_materials(components: dict[str, float]) -> dict[str, float]:
	""" take a dictionary denoting the molecular fractions of different materials in a mix, find any that are names of
	    gasses rather than actual nuclides, and expand those so that you have just the nuclides instead.  modify it in
	    place, but also return it
	"""
	initial_components = list(components.keys())
	for component in initial_components:
		try:
			parse_isotope_symbol(component)
		# if it's not a nuclide
		except ValueError:
			# then it must be a known compound material
			if component in LILAC_GAS_MATERIALS:
				# break it up into its constituent parts
				subcomponents = LILAC_GAS_MATERIALS[component].components
				for subcomponent in subcomponents.keys():
					components[subcomponent] = components.get(subcomponent, 0) + components[component]*subcomponents[subcomponent]
				del components[component]
			# if it's not known, give up
			else:
				raise ValueError(f"I don't know the material '{component}'")
	return components


def find_best_D3He_material_code(f3He: float) -> tuple[int, float]:
	""" find the LILAC material code that represents the mixture of D and ³He that most accurately
	    represents the specified fill ratio.  specifically, find the LILAC material with the nearest
	    atomic ³He fraction.
	    :return: the selected material code and the atomic ³He fraction that corresponds to it
	"""
	if 0 <= f3He <= 1:
		options = array(list(LILAC_D3He_MIXTURES.keys()))
		fractions = array([LILAC_D3He_MIXTURES[code] for code in options])
		best_index = argmin(abs(fractions - f3He))
		return options[best_index], fractions[best_index]
	raise ValueError(f"something's wrong with the ³He fill fraction {f3He:.0%}")


def get_solid_material_from_name(name: str) -> Material:
	""" create the Material object that describes the solid named by name.  see the [LILAC material
	    table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for
	    more information.
	"""
	# start by simply looking the name up in the dict of preset materials
	if name in LILAC_SOLID_MATERIALS:
		return LILAC_SOLID_MATERIALS[name]
	else:
		# otherwise try interpreting it as a material code
		try:
			material_code = int(name)
		except ValueError:
			raise IndexError(f"I don't recognize the material '{name}'")
		else:
			# search the preset materials for a matching one
			for material in LILAC_SOLID_MATERIALS.values():
				if material.material_code == material_code:
					return material
			# if there is no preset material, initialize one with the default options
			if material_code <= 92:
				return Material(material_code, eos=6, ionization=2, opacity=1)
			else:
				return Material(material_code, eos=8, ionization=1, opacity=1)


def get_gas_material_from_components(partial_pressures: dict[str, float]) -> Material:
	""" create the Material object that describes the solid named by name.  see the [LILAC material
	    table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for
	    more information.
	    :param partial_pressures: a dictionary where each key is the name of a nuclide and each
	                              value is its partial pressure in atm
	"""
	total_pressure = sum(partial_pressures.values())
	components = {symbol for symbol, pressure in partial_pressures.items() if pressure > 0}

	# convert the molecular partial pressures to atomic fractions
	atomic_fractions = defaultdict(lambda: 0.)
	for component in components:
		atomic_fractions[component] = partial_pressures[component]
		try:
			if parse_isotope_symbol(component)[0] in [1, 7, 8, 9, 17]:  # watch out for diatomic elements
				atomic_fractions[component] *= 2
		except ValueError:
			pass
	normalization = sum(atomic_fractions.values())
	for component in components:
		atomic_fractions[component] /= normalization

	for named_gas in LILAC_GAS_MATERIALS.keys():
		# if this is a particular gas that we have in the table
		if atomic_fractions[named_gas] > .9:
			# we can use that bilt-in material
			material = copy(LILAC_GAS_MATERIALS[named_gas])
			material.pressure = total_pressure
			if len(components) > 1:
				material.material_code *= -1  # marking it with a negative material code if it needs to be mixed
				material.components = drop_zeros(expand_compound_materials(atomic_fractions))
			return material

	atomic_fractions = expand_compound_materials(atomic_fractions)

	if sum(atomic_fractions[nuclide] for nuclide in ["¹H", "D", "T", "³He"]) > .9:
		# if this is a mixture of hydrogen isotopes
		if atomic_fractions["³He"] == 0:
			# we can use bilt-in gas #101 and control ratios using ptrit
			material = copy(LILAC_GAS_MATERIALS["D"])
			material.pressure = total_pressure
			material.components = drop_zeros(atomic_fractions)
			# if it has impurities (other than ³He), mark it with a negative material code
			if not components.issubset({"¹H", "D", "T"}):
				material.material_code *= -1
			return material

		# if this is a mixture of hydrogen isotopes and ³He
		else:
			# we can use a bilt-in gas iff there's a matcode that exactly matches the desired ratio
			best_material_code, found_3He_fraction = find_best_D3He_material_code(
				atomic_fractions["³He"])
			fraction_error = abs(found_3He_fraction - atomic_fractions["³He"])
			material = copy(LILAC_GAS_MATERIALS["D³He"])
			material.material_code = best_material_code
			material.pressure = total_pressure
			material.components = drop_zeros(atomic_fractions)
			if not components.issubset({"¹H", "D", "T", "³He"}) or fraction_error > 1e-2:
				# otherwise, mark this material with a negative matcode
				material.material_code = -best_material_code
			return material

	# if it's something else, we can't yet do it
	raise NotImplementedError("I don't autodetect materials by species unless it's mostly hydrogen and ³He.")
