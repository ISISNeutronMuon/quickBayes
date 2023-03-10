import unittest
import numpy as np
from quasielasticbayes.v2.fitting.gofit_engine import GoFitEngine
from quasielasticbayes.test_helpers.template_fit_test import FitEngineTemplate


class GoFitEngineTest(FitEngineTemplate, unittest.TestCase):

    @staticmethod
    def get_test_engine(x, y, e):
        engine = GoFitEngine(x, y, e,
                             lower=[-10, -10],
                             upper=[10, 10],
                             samples=10)
        engine.set_N_params(2)
        return engine

    @staticmethod
    def get_name():
        return "gofit"

    @staticmethod
    def get_basic_fit_params():
        return [0.968, 0.130], [0.046, 0.086]

    @staticmethod
    def get_chi_squared():
        return 2.497

    @staticmethod
    def get_covariance():
        return np.array([np.array([0.010, -0.005]),
                         np.array([-0.005, 0.003])])

    @staticmethod
    def get_basic_fit_values():
        expected_y = [.130, 1.098, 2.066, 3.034]
        expected_e = [0.113, 0.074, 0.074, 0.113]
        expected_d = [0.030, -0.102, 0.166, -0.116]
        expected_de = [0.151, 0.117, 0.133, 0.151]
        return expected_y, expected_e, expected_d, expected_de

    @staticmethod
    def get_spline_params():
        return [0.864, 0.100], [0.037, 0.022]

    @staticmethod
    def get_spline_fits():
        expected_y = [.100, 0.191, 0.282, 0.373, 0.464, 0.555,
                      0.646, 0.737, 0.828, 0.919]
        expected_e = [0.022, 0.019, 0.016, 0.013, 0.012, 0.012,
                      0.013, 0.015, 0.017, 0.020]
        expected_d = [-0.162, 0.049, -0.094, -0.185, -0.047, -0.165,
                      0.010, -0.139, -0.012, -0.038]
        expected_de = [0.055, 0.053, 0.052, 0.052, 0.051, 0.051,
                       0.052, 0.052, 0.053, 0.054]
        return expected_y, expected_e, expected_d, expected_de

    @staticmethod
    def get_low_stat_params():
        return [0.871, 0.174], [0.052, 0.029]

    @staticmethod
    def get_low_stat_fits():
        expected_y = [.174, 0.265, 0.357, 0.449, 0.540, 0.632,
                      0.724, 0.815, 0.907, 0.999]
        expected_e = [0.031, 0.027, 0.022, 0.019, 0.017, 0.017,
                      0.019, 0.022, 0.027, 0.031]
        expected_d = [-0.089, 0.123, -0.019, -0.110, 0.030, -0.088,
                      0.088, -0.061, 0.066, 0.042]
        expected_de = [0.059, 0.057, 0.055, 0.053, 0.053, 0.053,
                       0.053, 0.055, 0.057, 0.059]

        return expected_y, expected_e, expected_d, expected_de

    @staticmethod
    def get_spline_chi2():
        return {'low': 3.086, 'high': 5.389}

    @staticmethod
    def get_spline_covar():
        high = np.array([np.array([0.001, -0.001]),
                         np.array([-0.001, 0.0004])])
        low = np.array([np.array([0.003, -0.001]),
                        np.array([-0.001, 0.0004])])

        return {'high': high, 'low': low}

    """
     extra tests for scipy fit engine
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
    """


if __name__ == '__main__':
    unittest.main()
