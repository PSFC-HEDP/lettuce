import numpy as np
import numpy.testing as npt
import pytest
from numpy import array, concatenate, inf

from python.utilities import gradient, width, apparent_brightness, degrade_laser_pulse, to_superscript


def test_degrade_laser_pulse():
	npt.assert_equal(degrade_laser_pulse(array([0., 1., 1., 2., 0.]), 0.0),
	                 array([0., 1., 1., 2., 0.]))
	npt.assert_equal(degrade_laser_pulse(array([0., 1., 1., 2., 0.]), 0.5),
	                 array([0., 1., 1., 0., 0.]))
	npt.assert_equal(degrade_laser_pulse(array([0., 1., 1., 2., 0.]), 1.0),
	                 array([0., 0., 0., 0., 0.]))


def test_to_superscript():
	assert to_superscript("321") == "³²¹"
	with pytest.raises(KeyError):
		to_superscript("a")


def test_gradient():
	x = array([0., 1., 2.7, 3.6])
	y = array([1., 1., -0., 3.])
	np_gradient_dydx = np.gradient(y, x)
	npt.assert_equal(gradient(y, x), np_gradient_dydx)
	x = array([-2., 0., 1., 2.7, 3.6])
	y = array([1., 1., 1., 0., 3.])
	np_gradient_dydx = np.gradient(y, x)
	npt.assert_equal(gradient(y, x), concatenate([[0., 0.], np_gradient_dydx[2:]]))
	x = array([-2., 0., 1., 2.7, 3.6])
	y = array([3., 0., 1., 1., 1.])
	np_gradient_dydx = np.gradient(y, x)
	npt.assert_equal(gradient(y, x), concatenate([np_gradient_dydx[:-2], [0., 0.]]))


def test_width():
	npt.assert_allclose(width(array([0., 1., 2., 3.]),
	                          array([0., 1., 2., 0.])), 1.5)
	npt.assert_allclose(width(array([0., 1., 2., 3., 4.]),
	                          array([0., 2., 3., 1., 2.])), 2.0)
	npt.assert_allclose(width(array([0., 1., 2., 3., 4., 5.]),
	                          array([0., 2., 0., 2., 2., 0.])), 1.0)
	npt.assert_allclose(width(array([0., 1., 2., 3.]),
	                          array([2., 1., 0., 0.])), inf)
	npt.assert_allclose(width(array([0., 1., 2., 3.]),
	                          array([0., 0., 1., 2.])), inf)
	npt.assert_allclose(width(array([0., 1., 2., 3.]),
	                          array([0., 0., 0., 0.])), inf)
	npt.assert_allclose(width(array([0.]),
	                          array([1.])), inf)
	npt.assert_allclose(width(array([]),
	                          array([])), inf)


def test_apparent_brightness():
	Z = array(1.5)
	ne = array(1e+20)
	Te = array(3)
	baseline = apparent_brightness(Z, ne, Te, 0)
	dense = apparent_brightness(Z, 2*ne, Te, 0)
	assert dense > baseline
	empty = apparent_brightness(Z, 0*ne, Te, 0)
	assert empty == 0
	hot = apparent_brightness(Z, ne, 2*Te, 0)
	assert hot > baseline
	cold = apparent_brightness(Z, ne, 1e-20*Te, 0)
	assert cold == 0
	filtered = apparent_brightness(Z, ne, Te, 10)
	assert filtered < baseline
	with pytest.raises(ValueError):
		apparent_brightness(0*Z, ne, Te, 0)
	with pytest.raises(ValueError):
		apparent_brightness(Z, ne, 0*Te, 0)
