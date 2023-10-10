import numpy as np
import numpy.testing as npt
import pytest
from numpy import nan

from python.utilities import load_pulse_shape, find_best_D3He_material_code, get_shell_material_from_name, Material, \
	get_gas_material_from_components, parse_gas_components


def test_load_pulse_shape():
	time, power = load_pulse_shape("SG10v001", 7)
	npt.assert_allclose(np.sum(power*np.gradient(time)), 7)


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
	assert get_shell_material_from_name("DT") == Material(102, tritium_fraction=.50)
	assert get_shell_material_from_name("Be") == Material(4)
	with pytest.raises(IndexError):
		get_shell_material_from_name("surprise")


def test_get_gas_material_from_components():
	assert get_gas_material_from_components({"T": 1, "3He": 2}) == Material(295, tritium_fraction=.50, pressure=3)


def test_parse_gas_components():
	assert parse_gas_components("40atm Au +.3U") == {"Au": 40, "U": .3}
