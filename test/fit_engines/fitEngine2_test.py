import unittest
from numpy import ndarray
from typing import Callable
import numpy as np
from scipy.optimize import curve_fit
from quasielasticbayes.fitting.fit_engine import FitEngine
from quasielasticbayes.test_helpers.template_fit_test import FitEngineTemplate


class SimpleTestEngine(FitEngine):
    def __init__(self, x_data: ndarray, y_data: ndarray, e_data: ndarray):
        """
        Creates a test class using scipy
        :param x_data: the x data to fit
        :param y_data: the y data to fit
        :param e_data: the error data to fit
        """
        super().__init__("test", x_data, y_data, e_data)
        self._covar_cf = None

    def _do_fit(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                func: Callable) -> ndarray:
        """
        Calls scipy curve fit
        :param x_data: the x data to fit
        :param y_data: the y data to fit
        :param e_data: the error data to fit
        :return the fit parameters
        """
        params, covar = curve_fit(func, x_data, y_data,
                                  sigma=e_data, absolute_sigma=True)
        self._covar_cf = covar
        return params


class BaseFitEngineTest(FitEngineTemplate, unittest.TestCase):

    def get_test_engine(self, x, y, e):
        return SimpleTestEngine(x, y, e)

    def get_name(self):
        return "test"

    def get_basic_fit_params(self):
        return [0.986, 0.122], [0.047, 0.088]

    def get_chi_squared(self):
        return 2.347

    def get_covariance(self):
        return self.engine._covar_cf

    def get_basic_fit_values(self):
        expected_y = [.122, 1.108, 2.094, 3.081]
        expected_e = [0.116, 0.076, 0.076, 0.116]
        expected_d = [0.022, -0.092, 0.194, -0.069]
        expected_de = [0.153, 0.118, 0.134, 0.153]
        return expected_y, expected_e, expected_d, expected_de

    def get_spline_params(self):
        return [0.884, 0.095], [0.037, 0.022]

    def get_spline_fits(self):
        expected_y = [.095, 0.188, 0.281, 0.374, 0.467, 0.560,
                      0.653, 0.746, 0.839, 0.932]
        expected_e = [0.022, 0.019, 0.016, 0.013, 0.012, 0.012,
                      0.013, 0.015, 0.017, 0.020]
        expected_d = [-0.168, 0.046, -0.095, -0.185, -0.044, -0.160,
                      0.017, -0.130, -0.001, -0.025]
        expected_de = [0.055, 0.053, 0.052, 0.052, 0.051, 0.051,
                       0.052, 0.052, 0.053, 0.054]
        return expected_y, expected_e, expected_d, expected_de

    def get_low_stat_params(self):
        return [0.811, 0.204], [0.052, 0.029]

    def get_low_stat_fits(self):
        expected_y = [.204, 0.289, 0.375, 0.460, 0.545, 0.631,
                      0.716, 0.801, 0.887, 0.972]
        expected_e = [0.031, 0.027, 0.022, 0.019, 0.017, 0.017,
                      0.019, 0.022, 0.027, 0.031]
        expected_d = [-0.059, 0.147, -0.001, -0.099, 0.034, -0.089,
                      0.080, -0.075, 0.046, 0.015]
        expected_de = [0.059, 0.057, 0.055, 0.053, 0.053, 0.053,
                       0.053, 0.055, 0.057, 0.059]

        return expected_y, expected_e, expected_d, expected_de

    def get_spline_chi2(self):
        return {'low': 2.921, 'high': 5.366}

    def get_spline_covar(self):
        high = np.array([np.array([0.001, -0.001]),
                         np.array([-0.001, 0.0004])])
        low = np.array([np.array([0.003, -0.001]),
                        np.array([-0.001, 0.0004])])

        return {'high': high, 'low': low}


if __name__ == '__main__':
    unittest.main()
