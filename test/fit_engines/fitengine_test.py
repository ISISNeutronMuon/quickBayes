import unittest
from numpy import ndarray
from typing import Callable
import numpy as np
from scipy.optimize import curve_fit
from quasielasticbayes.v2.fitting.fit_engine import FitEngine
from quasielasticbayes.v2.functions.BG import LinearBG


def func(x_data):
    return 0.1 + 0.9*x_data


def data_1():
    return (np.array([0, 1, 2, 3]),
            np.array([0.1, 1.2, 1.9, 3.15]),
            np.array([0.1, 0.09, 0.11, 0.1]))


def data_2():
    x_data = np.linspace(0, 1, 20)
    np.random.seed(1)
    noise_stdev = 0.1
    y_data = np.random.normal(func(x_data), noise_stdev)
    e_data = 0.5*noise_stdev*np.ones(x_data.shape)
    return x_data, y_data, e_data


def data_3(x_data, y_data, e_data):
    xx = []
    yy = []
    ee = []

    for k in range(0, len(x_data), 2):
        xx.append(x_data[k])
        yy.append(y_data[k])
        ee.append(e_data[k])

    return np.array(xx), np.array(yy), np.array(ee)


# refactor checks and the above to a file in utils


class TestEngine(FitEngine):
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


class FitEngineTest(unittest.TestCase):

    def test_name(self):
        x_data, y_data, e_data = data_1()
        engine = TestEngine(x_data, y_data, e_data)
        self.assertEqual(engine.name, 'test')

    def test_fit_params(self):
        x_data, y_data, e_data = data_1()
        bg = LinearBG()

        engine = TestEngine(x_data, y_data, e_data)
        engine.do_fit(x_data, y_data, e_data, bg)

        params, errors = engine.get_fit_parameters()
        expected_p = [0.986, 0.122]
        expected_e = [0.047, 0.088]

        self.assertEqual(len(params), len(errors))
        for k in range(len(params)):
            self.assertAlmostEqual(params[k], expected_p[k], 3)
            self.assertAlmostEqual(errors[k], expected_e[k], 3)

    def test_fit_values(self):
        x_data, y_data, e_data = data_1()
        bg = LinearBG()

        engine = TestEngine(x_data, y_data, e_data)
        engine.do_fit(x_data, y_data, e_data, bg)

        xf, yf, ef, df, de = engine.get_fit_values()
        expected_y = [.122, 1.108, 2.094, 3.081]
        expected_e = [0.116, 0.076, 0.076, 0.116]
        expected_d = [0.022, -0.092, 0.194, -0.069]
        expected_de = [0.153, 0.118, 0.134, 0.153]

        self.assertEqual(len(xf), len(x_data))
        self.assertEqual(len(yf), len(x_data))
        self.assertEqual(len(ef), len(x_data))
        self.assertEqual(len(df), len(x_data))
        self.assertEqual(len(de), len(x_data))
        for k in range(len(xf)):
            self.assertAlmostEqual(xf[k], x_data[k], 3)
            self.assertAlmostEqual(yf[k], expected_y[k], 3)
            self.assertAlmostEqual(ef[k], expected_e[k], 3)
            self.assertAlmostEqual(df[k], expected_d[k], 3)
            self.assertAlmostEqual(de[k], expected_de[k], 3)

    def test_chi_squared(self):
        x_data, y_data, e_data = data_1()
        bg = LinearBG()

        engine = TestEngine(x_data, y_data, e_data)
        engine.do_fit(x_data, y_data, e_data, bg)
        self.assertAlmostEqual(engine.get_chi_squared(), 2.347, 3)

    def test_cov(self):
        # need to do more data to get an accurate covariance matrix
        x_data = np.linspace(0, 1, 10)
        y_data = func(x_data)
        e_data = 0.1*np.ones(len(x_data))
        bg = LinearBG()

        engine = TestEngine(x_data, y_data, e_data)
        engine.do_fit(x_data, y_data, e_data, bg)
        calculated = engine.get_covariance_matrix()
        scipy = engine._covar_cf

        self.assertEqual(len(calculated), len(scipy))
        self.assertEqual(len(calculated[0]), len(scipy[0]))

        for i in range(len(calculated)):
            for j in range(len(calculated[i])):
                self.assertAlmostEqual(calculated[i][j], scipy[i][j], 3)

    def test_spline_data_params(self):
        x_data, y_data, e_data = data_2()
        bg = LinearBG()
        xx, yy, ee = data_3(x_data, y_data, e_data)

        engine = TestEngine(xx, yy, ee)

        # fit with less data
        engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        engine.do_fit(x_data, y_data, e_data, bg)

        # check latest results
        params, errors = engine.get_fit_parameters()
        expected_p = [0.884, 0.095]
        expected_e = [0.037, 0.022]

        self.assertEqual(len(params), len(errors))
        for k in range(len(params)):
            self.assertAlmostEqual(params[k], expected_p[k], 3)
            self.assertAlmostEqual(errors[k], expected_e[k], 3)

        # check first results -> lower stats
        params, errors = engine.get_fit_parameters(0)
        expected_p = [0.811, 0.204]
        expected_e = [0.052, 0.029]

        self.assertEqual(len(params), len(errors))
        for k in range(len(params)):
            self.assertAlmostEqual(params[k], expected_p[k], 3)
            self.assertAlmostEqual(errors[k], expected_e[k], 3)

    def test_spline_data_fits(self):
        x_data, y_data, e_data = data_2()
        bg = LinearBG()
        xx, yy, ee = data_3(x_data, y_data, e_data)

        engine = TestEngine(xx, yy, ee)

        # fit with less data
        engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        engine.do_fit(x_data, y_data, e_data, bg)

        # check latest results
        xf, yf, ef, df, de = engine.get_fit_values()

        expected_y = [.095, 0.188, 0.281, 0.374, 0.467, 0.560,
                      0.653, 0.746, 0.839, 0.932]
        expected_e = [0.022, 0.019, 0.016, 0.013, 0.012, 0.012,
                      0.013, 0.015, 0.017, 0.020]
        expected_d = [-0.168, 0.046, -0.095, -0.185, -0.044, -0.160,
                      0.017, -0.130, -0.001, -0.025]
        expected_de = [0.055, 0.053, 0.052, 0.052, 0.051, 0.051,
                       0.052, 0.052, 0.053, 0.054]

        self.assertEqual(len(xf), len(xx))
        self.assertEqual(len(yf), len(xx))
        self.assertEqual(len(ef), len(xx))
        self.assertEqual(len(df), len(xx))
        self.assertEqual(len(de), len(xx))
        for k in range(len(xf)):
            self.assertAlmostEqual(xf[k], xx[k], 3)
            self.assertAlmostEqual(yf[k], expected_y[k], 3)
            self.assertAlmostEqual(ef[k], expected_e[k], 3)
            self.assertAlmostEqual(df[k], expected_d[k], 3)
            self.assertAlmostEqual(de[k], expected_de[k], 3)

        # check first results -> lower stats
        xf, yf, ef, df, de = engine.get_fit_values(0)

        expected_y = [.204, 0.289, 0.375, 0.460, 0.545, 0.631,
                      0.716, 0.801, 0.887, 0.972]
        expected_e = [0.031, 0.027, 0.022, 0.019, 0.017, 0.017,
                      0.019, 0.022, 0.027, 0.031]
        expected_d = [-0.059, 0.147, -0.001, -0.099, 0.034, -0.089,
                      0.080, -0.075, 0.046, 0.015]
        expected_de = [0.059, 0.057, 0.055, 0.053, 0.053, 0.053,
                       0.053, 0.055, 0.057, 0.059]

        self.assertEqual(len(xf), len(xx))
        self.assertEqual(len(yf), len(xx))
        self.assertEqual(len(ef), len(xx))
        self.assertEqual(len(df), len(xx))
        self.assertEqual(len(de), len(xx))
        for k in range(len(xf)):
            self.assertAlmostEqual(xf[k], xx[k], 3)
            self.assertAlmostEqual(yf[k], expected_y[k], 3)
            self.assertAlmostEqual(ef[k], expected_e[k], 3)
            self.assertAlmostEqual(df[k], expected_d[k], 3)
            self.assertAlmostEqual(de[k], expected_de[k], 3)

    def test_spline_chi_squared(self):
        x_data, y_data, e_data = data_2()
        bg = LinearBG()
        xx, yy, ee = data_3(x_data, y_data, e_data)

        engine = TestEngine(xx, yy, ee)

        # fit with less data
        engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        engine.do_fit(x_data, y_data, e_data, bg)

        self.assertAlmostEqual(engine.get_chi_squared(), 5.366, 3)
        self.assertAlmostEqual(engine.get_chi_squared(0), 2.921, 3)

    def test_spline_cov(self):
        x_data, y_data, e_data = data_2()
        bg = LinearBG()
        xx, yy, ee = data_3(x_data, y_data, e_data)

        engine = TestEngine(xx, yy, ee)

        # fit with less data
        engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        engine.do_fit(x_data, y_data, e_data, bg)

        calculated = engine.get_covariance_matrix()
        scipy = np.array([np.array([0.001, -0.001]),
                          np.array([-0.001, 0.0004])])

        self.assertEqual(len(calculated), len(scipy))
        self.assertEqual(len(calculated[0]), len(scipy[0]))

        for i in range(len(calculated)):
            for j in range(len(calculated[i])):
                self.assertAlmostEqual(calculated[i][j], scipy[i][j], 3)

        # check 1st fit
        calculated = engine.get_covariance_matrix(0)
        scipy = np.array([np.array([0.003, -0.001]),
                          np.array([-0.001, 0.0004])])

        self.assertEqual(len(calculated), len(scipy))
        self.assertEqual(len(calculated[0]), len(scipy[0]))

        for i in range(len(calculated)):
            for j in range(len(calculated[i])):
                self.assertAlmostEqual(calculated[i][j], scipy[i][j], 3)


if __name__ == '__main__':
    unittest.main()
