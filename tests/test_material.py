import pytest
from numpy import nan

from material import find_best_D3He_material_code, get_solid_material_from_name, Material, \
	get_gas_material_from_components, isotope_symbol, nuclide_symbol, parse_isotope_symbol, expand_compound_materials


def test_nuclide_symbol():
	assert nuclide_symbol(1, 1) == "p"
	assert nuclide_symbol(1, 2.014) == "d"
	assert nuclide_symbol(1, 3) == "t"
	assert nuclide_symbol(6, 12) == "¹²C"
	assert nuclide_symbol(6, 12.011) == "¹²C"


def test_isotope_symbol():
	assert isotope_symbol(1, 1) == "¹H"
	assert isotope_symbol(1, 2.014) == "D"
	assert isotope_symbol(1, 3) == "T"
	assert isotope_symbol(6, 12) == "¹²C"
	assert isotope_symbol(6, 12.011) == "natural C"


def test_parse_isotope_symbol():
	assert parse_isotope_symbol("D") == (1, 2.014)
	assert parse_isotope_symbol("¹²C") == (6, 12)
	assert parse_isotope_symbol("natural C") == (6, 12.011)
	assert parse_isotope_symbol("14.003C") == (6, 14.003)
	assert parse_isotope_symbol("C") == (6, 12.011)
	assert parse_isotope_symbol("232Th") == (90, 232)


def test_expand_compound_materials():
	assert expand_compound_materials({"D3He": .25, "DT": .25, "T": .5}) == {"D": .25, "T": .625, "³He": .125}
	with pytest.raises(ValueError):
		assert expand_compound_materials({"poop": 1})


def test_find_best_D3He_code():
	assert find_best_D3He_material_code(.00001) == (101, .0)
	assert find_best_D3He_material_code(.22) == (284, .2)
	assert find_best_D3He_material_code(.27) == (281, .28)
	assert find_best_D3He_material_code(.50) == (295, .5)
	assert find_best_D3He_material_code(.99999) == (103, 1.0)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(7)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(-1)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(nan)


def test_get_shell_material_from_name():
	assert get_solid_material_from_name("DT") == \
	       Material(102, eos=8, opacity=8, ionization=1, components={"T": 1/2})
	assert get_solid_material_from_name("Be") == \
	       Material(4, eos=6, opacity=1, ionization=2, components={})
	with pytest.raises(IndexError):
		get_solid_material_from_name("surprise")


def test_get_gas_material_from_components():
	assert get_gas_material_from_components({"T": 1, "³He": 2}) == \
	       Material(+295, eos=6, opacity=1, ionization=1, pressure=3, components={"T": 1/2, "³He": 1/2})
	assert get_gas_material_from_components({"T": 1, "³He": 2, "Ar": .04}) == \
	       Material(-295, eos=6, opacity=1, ionization=1, pressure=3.04, components={"T": 50/101, "³He": 50/101, "Ar": 1/101})
	assert get_gas_material_from_components({"T": 7, "³He": 2}) == \
	       Material(-286, eos=6, opacity=1, ionization=1, pressure=9, components={"T": 7/8, "³He": 1/8})
	assert get_gas_material_from_components({"D": 20.}) == \
	       Material(+101, eos=8, opacity=8, ionization=1, pressure=20., components={"D": 1})
	assert get_gas_material_from_components({"DT": 1}) == \
	       Material(+101, eos=8, opacity=8, ionization=1, pressure=1, components={"D": 1/2, "T": 1/2})
	assert get_gas_material_from_components({"air": 99, "Hg": 1}) == \
	       Material(-300, eos=4, opacity=1, ionization=1, pressure=100, components={"N": .99*.78, "O": .99*.21, "Ar": .99*.01, "Hg": .01})
