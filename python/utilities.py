import numpy as np
from matplotlib import pyplot as plt
from numpy import expand_dims, geomspace, arange, interp, exp, sqrt, isfinite, inf
from numpy.typing import NDArray
from scipy import integrate

from python import image_plate


def gradient(y: NDArray[float], x: NDArray[float], **kwargs):
	""" it's numpy.gradient but it fixes the roundoff errors you get with inequal steps """
	# find indices where there is exactly zero change
	flat = (y[0:-2] == y[1:-1]) & (y[1:-1] == y[2:])
	# account for the edges, conservatively assuming we plan to use 2nd order edges
	flat = np.concatenate([[flat[0]], flat, [flat[-1]]])
	# call the numpy function, slotting in zero where we know it should be zero
	return np.where(~flat, np.gradient(y, x, **kwargs), 0)


def width(x: NDArray[float], y: NDArray[float]):
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


def apparent_brightness(ionization: NDArray[float], electron_number_density: NDArray[float],
                        electron_temperature: NDArray[float], filter_stack=None, show_plot=False
                        ) -> NDArray[float]:
	""" how much of the emission would be detected by an image plate
		:param ionization: the spatio-temporally resolved average ion charge
		:param electron_number_density: the spacio-temporally resolved electron number density (cm^-3)
		:param electron_temperature: the spacio-temporally resolved electron temperature (keV)
		:param filter_stack: the filter specifications as (thickness (μm), material name)
		:param show_plot: whether to plot and show a spectrum before returning
		:return: the spacio-temporally resolved
	"""
	if filter_stack is None:
		filter_stack = []

	# account for sensitivity and transmission
	hν = expand_dims(geomspace(1e0, 1e3, 61), axis=tuple(1 + arange(ionization.ndim)))  # (keV)
	log_sensitivity = image_plate.log_xray_sensitivity(hν, filter_stack)

	# catch arithmetic errors before they happen
	if not np.all(isfinite(ionization) & (ionization > 0)):
		raise ValueError(f"some inputs to apparent_brightness() were invalid: Z={ionization:.5g}")
	if not np.all(isfinite(electron_number_density) & (electron_number_density >= 0)):
		raise ValueError(f"some inputs to apparent_brightness() were invalid: ne={electron_number_density:.5g}cm^-3")
	if not np.all(isfinite(electron_temperature) & (electron_temperature > 0)):
		raise ValueError(f"some inputs to apparent_brightness() were invalid: Te={electron_temperature:.5g}keV")

	# finally, account for the original spectrum (thus expanding the array to 3d)
	Z = ionization
	ne = electron_number_density
	ni = ne/Z
	Te = electron_temperature
	log_emission = -hν/Te
	emission = Z*ni*ne/sqrt(Te)*exp(log_emission + log_sensitivity)

	# plot the curve for my benefit
	if show_plot:
		plt.plot(hν, emission)
		plt.xscale("log")
		plt.xlabel(f"Energy (keV)")
		plt.ylabel(f"PSL (?)")
		plt.show()

	# finally, integrate over energy
	return integrate.trapz(x=hν, y=emission, axis=0)
