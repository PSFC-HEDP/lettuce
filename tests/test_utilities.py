import numpy as np
import numpy.testing as npt
import pytest
from numpy import array, concatenate, inf, ones

from python.utilities import gradient, width, apparent_brightness, degrade_laser_pulse, to_superscript, \
	select_key_indices, rebin


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


def test_select_key_indices():
	for inpoot, expected_output in [
		(array([0., 1., 3., 1., 1., 4., 0.]), (array([0, 2, 4, 5, 7]), array([1, 2, 4, 5]))),
		(array([0., 0., 0., 0.]), (array([0, 1, 2, 3, 4]), array([0, 1, 2, 3]))),
		(array([1., 1.]), (array([0, 1, 2]), array([0, 1]))),
		(array([]), (array([0]), array([]))),
	]:
		actual_output = select_key_indices(inpoot, 4)
		npt.assert_equal(actual_output[0], expected_output[0])
		npt.assert_equal(actual_output[1], expected_output[1])


def test_rebin():
	npt.assert_allclose(
		rebin(
			array([[6.], [2.], [7.], [3.], [1.], [8.], [5.], [3.], [0.], [9.]]),
			array([0, 3, 5, 10]), axis=0, weights=ones((10, 1))),
		array([[5.], [2.], [5.]]))
	with pytest.raises(ValueError):
		rebin(
			array([[6.], [2.], [7.], [3.], [1.], [8.], [5.], [3.], [0.], [9.]]),
			array([3, 5, 11]), axis=0, weights=ones((10, 1)))


def test_apparent_brightness():
	Z = array(1.5)
	ne = array(1e+20)
	Te = array(3)
	baseline = apparent_brightness(ne, Te, 0)
	dense = apparent_brightness(2*ne, Te, 0)
	assert dense > baseline
	empty = apparent_brightness(0*ne, Te, 0)
	assert empty == 0
	hot = apparent_brightness(ne, 2*Te, 0)
	assert hot > baseline
	cold = apparent_brightness(ne, 1e-20*Te, 1e-16)
	assert cold == 0
	filtered = apparent_brightness(ne, Te, 10)
	assert filtered < baseline
	with pytest.raises(ValueError):
		apparent_brightness(ne, 0*Te, 0)
