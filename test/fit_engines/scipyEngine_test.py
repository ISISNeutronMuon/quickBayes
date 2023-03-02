import unittest
import numpy as np
from quasielasticbayes.v2.fitting.scipy_engine import ScipyFitEngine
from quasielasticbayes.test_helpers.template_fit_test import FitEngineTemplate


class ScipyFitEngineTest(FitEngineTemplate, unittest.TestCase):

    def get_test_engine(self, x, y, e):
        return ScipyFitEngine(x, y, e,
                              lower=[-10, -10],
                              upper=[10, 10],
                              guess=[0, 0])

    def get_name(self):
        return "scipy"

    def get_basic_fit_params(self):
        return [0.986, 0.122], [0.045, 0.082]

    def get_chi_squared(self):
        return 2.347

    def get_covariance(self):
        return np.array([np.array([0.010, -0.005]),
                         np.array([-0.005, 0.003])])

    def get_basic_fit_values(self):
        expected_y = [.122, 1.108, 2.094, 3.081]
        expected_e = [0.108, 0.071, 0.073, 0.113]
        expected_d = [0.022, -0.092, 0.194, -0.069]
        expected_de = [0.147, 0.114, 0.132, 0.151]
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

    # extra tests for scipy fit engine
    def test_change_guess_and_bounds(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.set_guess_and_bounds([1, 1],
                                         [0, 0],
                                         [2, 2])
        self.assertEqual(self.engine._guess, [1, 1])
        self.assertEqual(self.engine._lower, [0, 0])
        self.assertEqual(self.engine._upper, [2, 2])

    def test_bad_guess(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1],
                                             [0, 0],
                                             [2, 2])

    def test_bad_lower(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1, 1],
                                             [0, 0, 0],
                                             [2, 2])

    def test_bad_upper(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1],
                                             [0, 0],
                                             [2])

    def test_bad_guess_and_lower(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1],
                                             [0],
                                             [2, 2])

    def test_bad_guess_and_upper(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1, 1, 1],
                                             [0, 0],
                                             [2, 2, 2])

    def test_bad_lower_and_upper(self):
        # not going to do a fit so data can be empty
        x_data, y_data, e_data = [], [], []
        self.engine = self.get_test_engine(x_data, y_data, e_data)
        with self.assertRaises(ValueError):
            self.engine.set_guess_and_bounds([1, 1],
                                             [0],
                                             [2])


if __name__ == '__main__':
    unittest.main()
