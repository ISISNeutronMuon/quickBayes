import unittest
import numpy as np
from quickBayes.functions.BG import LinearBG
from quickBayes.fitting.fit_utils import (log10_hessian_det,
                                          chi_squared,
                                          param_errors,
                                          derivative,
                                          fit_errors,
                                          var,
                                          res)


class FitUtilsTest(unittest.TestCase):

    def test_log10HessDet(self):
        covar = np.array([np.array([1, -5]), np.array([2, 1])])
        result = log10_hessian_det(covar)
        self.assertAlmostEqual(result, -1.041, 3)

    def test_log10HessDetBadData(self):
        # this data will produce a negative arg for a log
        covar = np.array([np.array([-.08, -0.04]), np.array([-0.04, 0.9])])
        result = log10_hessian_det(covar)
        # should return a negative number instead of NAN
        self.assertAlmostEqual(result, -9.0, 3)

    def test_chi2(self):
        x = np.array([0, 1, 2, 3, 4])
        y = np.array([-0.9, 1.1, 2.05, 2.8, 2.9])
        e = np.array([0.1, 0.1, 0.1, 0.2, 0.2])
        fit = x - 1.

        params = [1., -1]

        chi_2 = chi_squared(x, y, e, fit, params)

        self.assertAlmostEqual(chi_2, 82.833, 3)

    def test_param_errors(self):

        covar = np.array([np.array([0.9, -.08]), np.array([-.08, 1.02])])
        errors = param_errors(covar)
        self.assertAlmostEqual(errors[0], 0.949, 3)
        self.assertAlmostEqual(errors[1], 1.010, 3)

    def test_detivative(self):
        x = np.linspace(0, 5)

        def func(x, m, c):
            return m*x + c

        params = [1, -2]
        result = derivative(x, params, func)

        self.assertEqual(len(result[0]), len(x))
        for k in range(len(x)):
            self.assertAlmostEqual(result[0][k], x[k], 3)

        self.assertEqual(len(result[1]), len(x))
        for k in range(len(x)):
            self.assertAlmostEqual(result[1][k], 1.0, 3)

    def test_fit_errors(self):
        x = np.array([0, 1, 2, 3])
        y = 2*x + .1
        params = [2, .1]

        covar = np.array([np.array([1., -.8]), np.array([-.8, 1])])
        df_by_dp = [x, 1]

        errors = fit_errors(x, params, y, covar, df_by_dp)
        result = [1.321, 0.835, 1.772, 3.012]

        self.assertEqual(len(errors), len(result))
        for k in range(len(result)):
            self.assertAlmostEqual(errors[k], result[k], 3)

    def test_var(self):
        x = np.array([0, 1, 2, 3])
        y = 2*x + .1
        params = [1.9, .1]
        bg = LinearBG()
        result = var(bg, x, y, params)
        self.assertAlmostEqual(result, .14, 3)

    def test_res(self):
        x = np.array([0, 1, 2, 3])
        y = 2*x + .1
        e = 0.1*np.ones(len(y))
        params = [1.9, .1]
        bg = LinearBG()
        result = res(bg, x, y, e, params)
        self.assertAlmostEqual(result, 14., 3)


if __name__ == '__main__':
    unittest.main()
