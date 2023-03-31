import unittest
import numpy as np
from quickBayes.fitting.scipy_engine import ScipyFitEngine
from quickBayes.test_helpers.template_scipy_fit import ScipyFitTemplate


class ScipyFitEngineTest(ScipyFitTemplate, unittest.TestCase):

    @staticmethod
    def get_test_engine(x, y, e):
        return ScipyFitEngine(x, y, e,
                              lower=[-10, -10],
                              upper=[10, 10],
                              guess=[0, 0])

    @staticmethod
    def get_name():
        return "scipy"

    @staticmethod
    def get_basic_fit_params():
        return [0.986, 0.122], [0.045, 0.082]

    @staticmethod
    def get_covariance():
        return np.array([np.array([0.010, -0.005]),
                         np.array([-0.005, 0.003])])

    @staticmethod
    def get_basic_fit_values():
        expected_y = [.122, 1.108, 2.094, 3.081]
        expected_e = [0.108, 0.071, 0.073, 0.113]
        expected_d = [0.022, -0.092, 0.194, -0.069]
        expected_de = [0.147, 0.114, 0.132, 0.151]
        return expected_y, expected_e, expected_d, expected_de

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
