from __future__ import annotations

from copy import copy
from collections import defaultdict
from typing import Optional, Union

from utilities import to_superscript


class Material:
	def __init__(
			self, material_code: int, *,
			eos: int, opacity: Union[int, str], ionization: int,
			protium_fraction: float = 0.0, tritium_fraction: float = 0.0,
			density: Optional[float] = None, pressure: Optional[float] = None):
		""" define a material by its composition and density
		    :param material_code: the LILAC material code (see LILAC user guide)
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
		self.eos = eos
		self.opacity = opacity
		self.ionization = ionization

	def __eq__(self, other: Material):
		return self.material_code == other.material_code and \
		       self.protium_fraction == other.protium_fraction and \
		       self.tritium_fraction == other.tritium_fraction and \
		       self.density == other.density and \
		       self.pressure == other.pressure and \
		       self.eos == other.eos and \
		       self.opacity == other.opacity and \
		       self.ionization == other.ionization

	def __str__(self):
		return f"Material({self.material_code!r}, eos={self.eos!r}, opacity={self.opacity!r}, " \
		       f"ionization={self.ionization!r}, protium_fraction={self.protium_fraction!r}, " \
		       f"tritium_fraction={self.tritium_fraction!r}, density={self.density!r}, " \
		       f"pressure={self.pressure!r})"


# list of atomic symbols
ATOMIC_SYMBOLS = [
	"n", "H", "He",
	"Li", "Be", "B", "C", "N", "O", "F", "Ne",
	"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
	"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
	"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
	"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
	"Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
	"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf"]


# how to specify various common solids in LILAC
# see the [LILAC material table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for more information.
LILAC_SOLID_MATERIALS = {
	"H":           Material(102, eos=8, ionization=1, opacity=8, protium_fraction=1.0),
	"DT":          Material(102, eos=8, ionization=1, opacity=8, tritium_fraction=0.5),
	"Li":          Material(3, eos=6, ionization=2, opacity=1),
	"Be":          Material(4, eos=6, ionization=2, opacity=1),
	"diamond":     Material(6, eos=6, ionization=2, opacity=1, density=3.50),
	"HDC":         Material(6, eos=6, ionization=2, opacity=1, density=3.50),
	"Al":          Material(13, eos=6, ionization=2, opacity=1),
	"Si":          Material(14, eos=6, ionization=2, opacity=1),
	"Fe":          Material(26, eos=6, ionization=2, opacity=1),
	"Ge":          Material(32, eos=6, ionization=2, opacity=1),
	"Ta":          Material(73, eos=6, ionization=2, opacity=1),
	"Au":          Material(79, eos=6, ionization=4, opacity=21),
	"U":           Material(92, eos=6, ionization=2, opacity=1),
	"CH":          Material(110, eos=8, ionization=1, opacity=8, density=1.03),
	"strong CD":   Material(113, eos=8, ionization=1, opacity=8, density=1.09),
	"CD":          Material(115, eos=8, ionization=1, opacity=1, density=1.10),
	"SiO2":        Material(150, eos=6, ionization=1, opacity=1, density=2.20),
	"glass":       Material(150, eos=6, ionization=1, opacity=1, density=2.20),
	"polystyrene": Material(510, eos=8, ionization=1, opacity=8, density=1.05),
}

# how to specify various common gases in LILAC
LILAC_GAS_MATERIALS = {
	"D": Material(101, eos=8, ionization=1, opacity=8),
	"DT": Material(101, eos=8, ionization=1, opacity=8, tritium_fraction=0.5),
	"He": Material(2, eos=6, ionization=2, opacity=1),
	"N": Material(7, eos=6, ionization=2, opacity=1),
	"Ne": Material(10, eos=6, ionization=2, opacity=1),
	"Ar": Material(18, eos=8, ionization=1, opacity=20),
	"D3He": Material(295, eos=6, ionization=1, opacity=1),
	"air": Material(300, eos=4, ionization=1, opacity=1),
}

# material codes that correspond to mixtures of D3He
LILAC_D3He_MIXTURES = {
	103: 1.00000, 289: 0.99000, 278: 0.92308, 292: 0.90000,
	280: 0.87500, 298: 0.80000, 291: 0.70000, 297: 0.66667,
	296: 0.60000, 295: 0.50000, 283: 0.42857, 288: 0.40000,
	290: 0.33333, 279: 0.28571, 281: 0.28000, 284: 0.20000,
	282: 0.15000, 286: 0.11111, 287: 0.10000, 285: 0.07692,
	277: 0.02500, 299: 0.00990, 101: 0.00000,
}


def nuclide_symbol(atomic_number: int, mass_number: int) -> str:
	""" succinctly describe this particular nuclide """
	if atomic_number == 1 and mass_number == 2:
		return "D"
	elif atomic_number == 1 and mass_number == 3:
		return "T"
	else:
		return to_superscript(str(mass_number)) + ATOMIC_SYMBOLS[atomic_number]


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


def get_solid_material_from_name(name: str) -> Material:
	""" create the Material object that describes the solid named by name.  see the [LILAC material
	    table](https://lle-prod-gitlab.lle.rochester.edu/lilac/lilac/-/wikis/material-table) for
	    more information.
	"""
	if name in LILAC_SOLID_MATERIALS:
		return LILAC_SOLID_MATERIALS[name]
	else:
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
	atomic_fractions = defaultdict(lambda: 0.)
	for nuclide in nuclides:
		if nuclide in ["H", "D", "T", "N", "O", "F", "Cl"]:
			atomic_fractions[nuclide] = 2*partial_pressures[nuclide]
		else:
			atomic_fractions[nuclide] = partial_pressures[nuclide]
	normalization = sum(atomic_fractions.values())
	for nuclide in nuclides:
		atomic_fractions[nuclide] /= normalization

	for named_gas in LILAC_GAS_MATERIALS.keys():
		# if this is a particular gas that we have in the table
		if atomic_fractions[named_gas] > .9:
			print(f"ignoring anything in the gas other than {named_gas}")
			# we can use that bilt-in material
			material = copy(LILAC_GAS_MATERIALS[named_gas])
			material.pressure = total_pressure
			return material

	# if this is just a mixture of hydrogen isotopes
	if sum(atomic_fractions[nuclide] for nuclide in ["H", "D", "T"]) > .9:
		print("ignoring anything in the gas other than H, D, and T")
		# we can use bilt-in gas #101 and control ratios using ptrit
		# TODO: if you ever implement mixing, material -101 should use opacity table "/lle/data/opacity_tables/HDT_plt_50x50Dt_101.prp"
		material = copy(LILAC_GAS_MATERIALS["D"])
		material.protium_fraction = atomic_fractions["H"]
		material.tritium_fraction = atomic_fractions["T"]
		material.pressure = total_pressure
		return material

	# if this is a mixture of hydrogen and helium-3
	if sum(atomic_fractions[nuclide] for nuclide in ["H", "D", "T", "3He"]) > .9:
		print("ignoring anything in the gas other than H, D, T and 3He")
		# we can use a bilt-in gas iff there's a matcode that exactly matches the desired ratio
		material = copy(LILAC_GAS_MATERIALS["D3He"])
		material.protium_fraction = atomic_fractions["H"]
		material.tritium_fraction = atomic_fractions["T"]
		material.pressure = total_pressure
		best_material_code, found_atomic_fraction = find_best_D3He_material_code(atomic_fractions["3He"])
		if abs(found_atomic_fraction - atomic_fractions["3He"]) < 1e-2:
			material.material_code = best_material_code
			return material
		else:
			# otherwise, we will need to specify this material manually using `mater`
			print("warning: I haven't implemented the mater namelist")
			material.material_code = best_material_code
			return material

	# if it's something else, we can't yet do it
	raise NotImplementedError("I don't autodetect materials by species unless it's a simple gas mix.")
