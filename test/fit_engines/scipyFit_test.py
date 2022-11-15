import unittest
from numpy import ndarray
import numpy as np
from quasielasticbayes.v2.fitting.scipy_fit import scipy_curve_fit
from quasielasticbayes.v2.functions.gaussian import Gaussian


def mock_data(x: ndarray, a: float, mu: float,
              sigma: float) -> (ndarray, ndarray):
    """
    Produces a gaussian plus noise for testing fits
    :param x: x values to evaluate function over
    :param a: amplitude of gaussian
    :param mu: mean of gaussin
    :param sigma: sigma value of gaussin
    :return gaussian with nosie and errors
    """
    g = Gaussian()
    y_data = g(x, a, mu, sigma)
    np.random.seed(10)  # make sure rand no are always the same
    e = np.random.rand(len(x))
    errors = 0.2*y_data
    return y_data*(1. + 0.2*(e - 0.5)), errors


class ScipyFitTest(unittest.TestCase):

    def test_gaussian_fit(self):
        x = np.linspace(-2., 2., 100)
        amp = 0.83
        mu = 0.39
        sigma = 0.2
        y, e = mock_data(x, amp, mu, sigma)
        g = Gaussian()
        guess = g.get_guess()
        bounds = g.get_bounds()

        (chi2, hess_det,
         params, fit) = scipy_curve_fit(x, y, e, g, guess, bounds[0],
                                        bounds[1])

        self.assertAlmostEqual(params[0], amp, 2)
        self.assertAlmostEqual(params[1], mu, 2)
        self.assertAlmostEqual(params[2], sigma, 2)

        expect = g(x, *params)
        for j in range(len(fit)):
            self.assertAlmostEqual(fit[j], expect[j], 3)

        self.assertAlmostEqual(hess_det, 17.643, 3)
        self.assertAlmostEqual(chi2, 0.078, 3)


if __name__ == '__main__':
    unittest.main()
