import numpy as np
import numpy.testing as npt
import pytest
from pandas import isnull

from data_io import load_pulse_shape, \
	parse_gas_components, load_inputs_table, load_outputs_table, log_message, \
	fill_in_template, write_row_to_outputs_table


def test_load_inputs_table():
	inputs_table = load_inputs_table()
	assert inputs_table.index.dtype == str
	assert inputs_table["flux limiter"].dtype == float


def test_load_outputs_table():
	outputs_table = load_outputs_table()
	assert outputs_table["yield"].dtype == float


def test_write_row_to_outputs_table():
	with pytest.raises(KeyError):
		write_row_to_outputs_table({"yield": 4e17})
	with pytest.raises(KeyError):
		write_row_to_outputs_table({"name": "test", "x": 0, "y": 2})
	write_row_to_outputs_table({"name": "test", "yield": 4e17}, drop_previous_data=True)
	assert load_outputs_table().loc["test", "yield"] == 4e17
	assert isnull(load_outputs_table().loc["test", "bang-time"])
	write_row_to_outputs_table({"name": "test", "yield": 8e17, "bang-time": 1.5})
	assert load_outputs_table().loc["test", "yield"] == 8e17
	assert load_outputs_table().loc["test", "bang-time"] == 1.5


def test_log_message():
	log_message("testing.")


def test_fill_in_template():
	expected_output = "krabs is a nice person\n" \
	                  "0: electrode\n" \
	                  "1: diglett\n" \
	                  "2: nidoran\n" \
	                  "3: mankey\n" \
	                  "my favorite is diglett\n"
	actual_output = fill_in_template(
		"../../tests/template.txt",
		parameters={"noun": "nice person", "pokemon": ["electrode", "diglett", "nidoran", "mankey"]},
		flags={"true": True, "false": False},
		loops={"i": range(4)})
	assert actual_output == expected_output


def test_load_pulse_shape():
	time, power = load_pulse_shape("SG10v001", 7)
	npt.assert_allclose(np.sum(power*np.gradient(time)), 7)


def test_parse_gas_components():
	assert parse_gas_components("40atm Au +.3 235U") == {"natural Au": 40, "²³⁵U": .3}
	assert parse_gas_components("40atm air") == {"air": 40}
	with pytest.raises(ValueError):
		parse_gas_components("10atm D + 10atm 2.014H")
