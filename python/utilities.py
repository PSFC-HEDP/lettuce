import json
from datetime import datetime
from re import fullmatch, sub, DOTALL, search
from typing import Optional

import numpy as np
from numpy import isfinite
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

	# start with the parameter values
	for key, value in parameters.items():
		if value == "nan":
			raise ValueError("you should never pass 'nan' into an input deck.  is this pandas's doing?  god, I "
			                 "wish pandas wouldn't use nan as a missing placeholder; that's so incredibly not "
			                 "what it's for.  onestly I just wish nan didn't exist.   it's like javascript null.  "
			                 "anyway, check your inputs.  make sure none of them are empty or nan (except the "
			                 "last three, which may be empty).")
		if search(f"<<{key}>>", content):
			content = sub(f"<<{key}>>", value, content)
		else:
			raise KeyError(f"the parameter <<{key}>> was not found in the template {template_filename}.")

	# then do the flags
	if flags is not None:
		for key, active in flags.items():
			# if it's active, remove the <<if>> and <<endif>> lines
			if active:
				content = sub(f"<<if {key}>>\n", "", content)
				content = sub(f"<<endif {key}>>\n", "", content)
			# if it's inactive, remove the <<if>> and <<endif>> lines and everything between them
			else:
				content = sub(f"<<if {key}>>.*<<endif {key}>>\n", "", content, flags=DOTALL)
	# check to make sure we got it all
	remaining_blank = search("<<.*>>", content)
	if remaining_blank:
		raise KeyError(f"you tried to fill out the template {template_filename} without specifying "
		               f"the value of <<{remaining_blank.group()}>>")

	return content


def load_pulse_shape(pulse_shape_name: str, total_energy: float) -> tuple[NDArray[float], NDArray[float]]:
	""" load a pulse shape from disk
	    :return: the time (ns) and total laser power (TW)
	"""
	if not isfinite(total_energy):
		raise ValueError("any arithmetic operation involving nan should raise an error in my opinion. since "
		                 "numpy won't do that for me, I'm doing it here now. this laser energy is nan. check "
		                 "your input table before you wreck your input table.")
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
