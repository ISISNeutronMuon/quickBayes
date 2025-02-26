import unittest
from numpy import ndarray
from typing import Callable
from scipy.optimize import curve_fit
from quickBayes.fitting.fit_engine import FitEngine
from quickBayes.test_helpers.template_scipy_fit import ScipyFitTemplate


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


class BaseFitEngineTest(ScipyFitTemplate, unittest.TestCase):

    @staticmethod
    def get_test_engine(x, y, e):
        return SimpleTestEngine(x, y, e)

    @staticmethod
    def get_name():
        return "test"

    @staticmethod
    def get_basic_fit_params():
        return [0.986, 0.122], [0.047, 0.088]

    def get_covariance(self):
        return self.engine._covar_cf

    @staticmethod
    def get_basic_fit_values():
        expected_y = [.122, 1.108, 2.094, 3.081]
        expected_e = [0.116, 0.076, 0.076, 0.116]
        expected_d = [0.022, -0.092, 0.194, -0.069]
        expected_de = [0.153, 0.118, 0.134, 0.153]
        return expected_y, expected_e, expected_d, expected_de


if __name__ == '__main__':
    unittest.main()
