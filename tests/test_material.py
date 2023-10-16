import pytest
from numpy import nan

from python.material import find_best_D3He_material_code, get_solid_material_from_name, Material, \
	get_gas_material_from_components


def test_find_best_D3He_code():
	assert find_best_D3He_material_code(.00001) == (101, .0)
	assert find_best_D3He_material_code(.22) == (284, .2)
	assert find_best_D3He_material_code(.27) == (284, .2)
	assert find_best_D3He_material_code(.50) == (295, .5)
	assert find_best_D3He_material_code(.99999) == (103, 1.0)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(7)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(-1)
	with pytest.raises(ValueError):
		find_best_D3He_material_code(nan)


def test_get_shell_material_from_name():
	assert get_solid_material_from_name("DT") == Material(102, eos=8, opacity=8, ionization=1, tritium_fraction=.50)
	assert get_solid_material_from_name("Be") == Material(4, eos=6, opacity=1, ionization=2)
	with pytest.raises(IndexError):
		get_solid_material_from_name("surprise")


def test_get_gas_material_from_components():
	assert get_gas_material_from_components({"T": 1, "3He": 2}) == Material(295, eos=6, opacity=1, ionization=1, tritium_fraction=.50, pressure=3)
	assert get_gas_material_from_components({"D": 1}) == Material(101, eos=8, opacity=8, ionization=1, pressure=1)
	assert get_gas_material_from_components({"air": 1, "cyanide": 0}) == Material(300, eos=4, opacity=1, ionization=1, pressure=1)
