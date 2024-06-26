from datetime import datetime
from typing import Union, Sequence, Dict, Tuple

import numpy as np
from fpdf import FPDF
from numpy import arange, interp, exp, sqrt, isfinite, inf, cumsum, where, concatenate, \
	linspace, unique, \
	zeros, empty, floor, histogram, float64, ndarray


def degrade_laser_pulse(original_pulse: ndarray, factor: float) -> ndarray:
	""" take a laser pulse and reduce its energy by the given factor by setting stuff at the end to zero. """
	time = arange(len(original_pulse))
	energy_integral = cumsum(original_pulse)/np.sum(original_pulse)
	cutoff_index = interp(1 - factor, energy_integral, time)
	return where(time <= cutoff_index, original_pulse, 0)


def to_superscript(string: str) -> str:
	""" convert a string to use the special unicode superscript letters """
	superscript = {
		"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
		"5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹",
	}
	return "".join(superscript.get(character, character) for character in string)


def from_superscript(string: str) -> str:
	""" convert a string to not use special unicode superscript letters """
	normalscript = {
		"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
		"⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
	}
	return "".join(normalscript.get(character, character) for character in string)


def drop_zeros(counts: Dict[str, Union[int, float]]) -> Dict[str, Union[int, float]]:
	""" copy a dictionary with all zero elements removed """
	output = {}
	for key in counts.keys():
		if counts[key] != 0:
			output[key] = counts[key]
	return output


def gradient(y: np.ndarray, x: np.ndarray, **kwargs):
	""" it's numpy.gradient but it fixes the roundoff errors you get with inequal steps """
	# find indices where there is exactly zero change
	flat = (y[0:-2] == y[1:-1]) & (y[1:-1] == y[2:])
	# account for the edges, conservatively assuming we plan to use 2nd order edges
	flat = np.concatenate([[flat[0]], flat, [flat[-1]]])
	# call the numpy function, slotting in zero where we know it should be zero
	return np.where(~flat, np.gradient(y, x, **kwargs), 0)


def width(x: np.ndarray, y: ndarray):
	""" calculate the full-width at half-max of a curve.
	    this function isn't very fast, but that's okay.
	    :return: the FWHM if calculable.  if the peak goes off the end or there is no peak, return inf.
	"""
	if x.shape != y.shape or x.ndim != 1:
		raise ValueError
	if x.size == 0:
		return inf
	maximum = np.argmax(y)
	height = y[maximum]
	if height <= 0:
		return inf
	x_left = None
	for i in range(maximum - 1, -1, -1):
		if y[i] <= height/2:
			x_left = interp(height/2, [y[i], y[i + 1]], [x[i], x[i + 1]])
			break
	x_rite = None
	for i in range(maximum + 1, x.size):
		if y[i] <= height/2:
			x_rite = interp(height/2, [y[i], y[i - 1]], [x[i], x[i - 1]])
			break
	if x_left is None or x_rite is None:
		return inf
	else:
		return x_rite - x_left


def select_key_indices(weights: ndarray, num_regions: int) -> Tuple[ndarray, ndarray]:
	""" given some 1D data representing a distribution, select a set of n indices for n points in that
	    distribution that characterize it well.  specifically, divide it into n regions of roughly equal
	    integral, then find the indices of the bin that contains the median of each region.
	    :return: num_regions+1 division indices, starting with 0 and ending with weights.size and having no duplicates, and
	             num_regions central indices, each falling within [division[i], division[i+1])
	"""
	if np.any(weights < 0):
		raise ValueError("negative inputs to this function cause ambiguity")
	if weights.size == 0:
		return zeros(1, dtype=int), zeros(0, dtype=int)

	x = arange(weights.size + 1)
	y = concatenate([[0], cumsum(weights)])
	# in the event of zero weights, act like the weights are uniform instead
	if np.min(y) == np.max(y):
		y = x

	divisions = interp(linspace(np.min(y), np.max(y), num_regions + 1), y, x)
	# snap them to integer values
	divisions = np.round(divisions).astype(int)
	# drop any duplicates (this means the size of the array won't always be num_regions+1 but that's fine)
	divisions = unique(divisions)
	# if the first few elements are 0, the first bin can be ambiguus.  snap it to 0.
	if divisions[0] != 0:
		assert y[divisions[0]] == 0
	divisions[0] = 0

	# finally, find the median point of each region
	centers = empty(divisions.size - 1, dtype=int)
	for i in range(divisions.size - 1):
		centers[i] = floor(interp((y[divisions[i]] + y[divisions[i + 1]])/2, y, x))

	return divisions, centers


def rebin(values: ndarray, bin_indices: Sequence[int], weights: ndarray, axis=0) -> ndarray:
	""" reduce the size of an array by averaging its values within certain index bins
	    :param values: the original values to be averaged with each other
	    :param bin_indices: the index at each bin edge. a value `values[i]` belongs to the bin `j` iff i ∈ [bin_indices[j], bin_indices[j + 1]).
	    :param weights: the weighting from each value
	    :param axis: the axis along which to apply the bins
	    :return: an array of averages, where the jth element is the average of all values in bin j
	"""
	if bin_indices[0] != 0:
		raise ValueError(f"the left edge of the first bin must be 0, not {bin_indices[0]}")
	elif bin_indices[-1] != values.shape[axis]:
		raise ValueError(f"the right edge of the last bin must be the size of values ({values.size}), not {bin_indices[-1]}")
	elif axis != 0 or values.ndim != 2:
		raise NotImplementedError("I haven't generalized this, sorry.")
	new_values = empty((len(bin_indices) - 1,) + values.shape[1:])
	for i in range(values.shape[1]):
		# use double precision because histogram is super sensitive to roundoff for some reason
		totals, _ = histogram(arange(values.shape[axis]), bins=bin_indices, weights=(weights*values)[:, i].astype(float64))
		numbers, _ = histogram(arange(values.shape[axis]), bins=bin_indices, weights=weights[:, i])
		new_values[:, i] = totals/numbers
	return new_values


def apparent_brightness(electron_number_density: ndarray,
                        electron_temperature: ndarray, energy_cutoff: float
                        ) -> ndarray:
	""" how much of the emission would be detected by an image plate
		:param electron_number_density: the spacio-temporally resolved electron number density (cm^-3)
		:param electron_temperature: the spacio-temporally resolved electron temperature (keV)
		:param energy_cutoff: x-ray energies will be integrated from this to positive infinity
		:return: the spacio-temporally resolved power per unit volume (units unknown)
	"""
	# catch arithmetic errors before they happen
	if not np.all(isfinite(electron_number_density) & (electron_number_density >= 0)):
		raise ValueError(f"some inputs to apparent_brightness() were invalid: ne = {np.min(electron_number_density):.5g} cm^-3")
	if not np.all(isfinite(electron_temperature) & (electron_temperature > 0)):
		raise ValueError(f"some inputs to apparent_brightness() were invalid: Te = {np.min(electron_temperature):.5g} keV")

	# finally, account for the original spectrum (thus expanding the array to 3d)
	return (electron_number_density**2 *
	        sqrt(electron_temperature) *
	        exp(-energy_cutoff/electron_temperature))


def create_pdf(title: str) -> FPDF:
	""" start an FPDF object a nicely formatted header """
	pdf = FPDF("landscape", "mm", "A4")
	pdf.add_font("Noto", "", "resources/fonts/NotoSans-Light.ttf", uni=True)
	pdf.add_font("Noto", "B", "resources/fonts/NotoSans-Bold.ttf", uni=True)
	pdf.add_page()
	pdf.set_font("Noto", "B", 20)
	pdf.set_title(title)
	pdf.write(15, title)
	# put the date in the upper corner
	pdf.set_font("Noto", "", 12)
	pdf.set_x(-55)
	pdf.write(15, datetime.today().strftime("%Y %b %d %H:%M:%S"))
	pdf.ln()
	return pdf



class RecordNotFoundError(Exception):
	pass

class InvalidSimulationError(Exception):
	pass
