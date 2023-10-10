import logging
from argparse import ArgumentParser

from matplotlib import pyplot as plt
from pandas import read_csv

from python.utilities import load_pulse_shape, get_shell_material_from_name, get_gas_material_from_components, \
	parse_gas_components


def start_lilac_run(name: str):
	""" arrange the inputs for a LILAC run in a machine-readable format and submit it to slurm
	    :param name: a string that will let us look up this run's inputs in the run_inputs.csv table
	"""
	inputs_table = read_csv("run_inputs.csv", skipinitialspace=True, index_col="name", dtype={"name": str})
	try:
		inputs = inputs_table.loc[name]
	except KeyError:
		logging.error(f"please add shot '{name}' to the run_inputs.csv file.")
		return
	if inputs.ndim != 1:
		logging.error(f"run '{name}' appears more than once in the input/shot_info.csv file!")
		return

	time, power = load_pulse_shape(inputs["pulse shape"], inputs["laser energy"])
	shell_material = get_shell_material_from_name(inputs["shell material"])
	fill_material = get_gas_material_from_components(parse_gas_components(inputs["fill"]))

	print(shell_material)
	print(fill_material)
	plt.figure()
	plt.plot(time, power)
	plt.show()

if __name__ == "__main__":
	parser = ArgumentParser(
		prog="start_lilac_run.sh",
		description = "Arrange the inputs for a LILAC run in a machine-readable format and submit it to slurm")
	parser.add_argument(
		"name", type=str,
		help="The name of the run, as specified in run_inputs.csv")
	args = parser.parse_args()

	start_lilac_run(args.name)
