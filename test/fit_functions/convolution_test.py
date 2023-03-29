import unittest
from numpy import ndarray
import numpy as np
from quickBayes.functions.gaussian import Gaussian
from quickBayes.functions.convolution import (
        ConvolutionWithResolution as conv)
from quickBayes.utils.crop_data import crop


def analytic(x: ndarray, amp: float, mu: float, sig: float,
             r_mu: float, r_sig: float) -> ndarray:
    """
    Expected results from wofram site on convolution
    """
    expect = np.exp(- pow(x-(mu + r_mu), 2)/(2.*(r_sig**2 + sig**2)))
    expect *= amp/np.sqrt(2.*np.pi*(r_sig**2 + sig**2))
    return expect


class ConvolutionTest(unittest.TestCase):

    def test_update_resolution(self):
        def func(x):
            return x*x - 3.*x + 1.2

        x = np.linspace(-5, 5, 6)
        y = func(x)
        c = conv(x, y, -6, 6)

        new_x = np.linspace(-5, 5, 100)
        c.update_x_range(new_x)
        ry = c._ry
        expect = func(new_x)
        # normalise the expected values:
        expect /= sum(expect)

        self.assertEqual(len(ry), len(new_x))
        for j in range(len(ry)):
            self.assertAlmostEqual(ry[j], expect[j], 3)

    def test_conv_call(self):
        """
        Need the x range to go to zero at the ends to
        prevent edge effects
        Also need lots of data points as
        the shape needs to be well defined
        """
        x = np.linspace(-15., 15, 9000)
        """
        need to set amp=1
        Since the kernel is normalised
        Let the other fit parameters
        absorb the value
        """
        res = Gaussian()
        r_a = 1.0
        r_mu = 4.0
        r_sig = 1.1
        res_y = res(x, r_a, r_mu, r_sig)

        c = conv(x, res_y, -16.0, 16.)

        g = Gaussian()
        c.add_function(g)

        amp = 1.
        mu = -2.4
        sig = 0.8

        # expected results from wofram site on convolution
        expect = analytic(x, amp, mu, sig, r_mu, r_sig)
        y = c(x, amp, mu, sig)

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_conv_call_with_crop(self):
        """
        Need the x range to go to zero at the ends to
        prevent edge effects
        Also need lots of data points as
        the shape needs to be well defined
        """
        x = np.linspace(-15., 15, 9000)
        """
        need to set amp=1
        Since the kernel is normalised
        Let the other fit parameters
        absorb the value
        """
        res = Gaussian()
        r_a = 1.0
        r_mu = 4.0
        r_sig = 1.1
        res_y = res(x, r_a, r_mu, r_sig)

        c = conv(x, res_y, -10.0, 10.)

        g = Gaussian()
        c.add_function(g)

        amp = 1.
        mu = -2.4
        sig = 0.8
        x, _, _ = crop(x, res_y, None, -10, 10)

        expect = analytic(x, amp, mu, sig, r_mu, r_sig)
        y = c(x, amp, mu, sig)

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_conv_report(self):
        report = {"old": [1]}
        x = np.linspace(0, 1)
        g = Gaussian()
        c = conv(x, x, 0, 1)
        c.add_function(g)
        out = c.report(report, 3.2, -1, 2.5)

        self.assertEqual(out["f1.f1.Amplitude"], [3.2])
        self.assertEqual(out["f1.f1.Mean"], [-1])
        self.assertEqual(out["f1.f1.Sigma"], [2.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 4)

    def test_conv_report2(self):
        report = {"old": [1]}
        x = np.linspace(0, 1)
        g = Gaussian()
        g2 = Gaussian()
        c = conv(x, x, 0, 1)
        c.add_function(g)
        c.add_function(g2)
        out = c.report(report, 3.2, -1, 2.5, 5, 6, 7)

        self.assertEqual(out["f1.f1.Amplitude"], [3.2])
        self.assertEqual(out["f1.f1.Mean"], [-1])
        self.assertEqual(out["f1.f1.Sigma"], [2.5])
        self.assertEqual(out["f1.f2.Amplitude"], [5])
        self.assertEqual(out["f1.f2.Mean"], [6])
        self.assertEqual(out["f1.f2.Sigma"], [7])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 7)

    def test_read(self):
        report = {"old": [1]}
        x = np.linspace(0, 1)
        g = Gaussian()
        c = conv(x, x, 0, 1)
        c.add_function(g)
        out = c.report(report, 3.2, -1, 2.5)
        params = c.read_from_report(out, 0)

        self.assertEqual(params, [3.2, -1, 2.5])

    def test_read2(self):
        report = {"old": [1]}
        x = np.linspace(0, 1)
        g = Gaussian()
        g2 = Gaussian()
        c = conv(x, x, 0, 1)
        c.add_function(g)
        c.add_function(g2)
        out = c.report(report, 3.2, -1, 2.5, 5, 6, 7)
        params = c.read_from_report(out, 0)
        self.assertEqual(params, [3.2, -1, 2.5, 5, 6, 7])

    def test_N_params(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        self.assertEqual(c.N_params, 0)
        g = Gaussian()
        c.add_function(g)
        self.assertEqual(c.N_params, 3)

    def test_guess(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        g = Gaussian()
        c.add_function(g)
        self.assertEqual(c.get_guess(), [1., 0., 0.1])

    def test_bounds(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        g = Gaussian()
        c.add_function(g)
        lower, upper = c.get_bounds()
        self.assertEqual(lower, [0., -1., 0.0])
        self.assertEqual(upper, [np.inf, 1., np.inf])

    def test_set_guess(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        g = Gaussian()
        g2 = Gaussian()
        c.add_function(g)
        c.add_function(g2)
        self.assertEqual(c.get_guess(), [1., 0., 0.1, 1., 0, 0.1])

        c.set_guess([2, 3, 4], 0)
        self.assertEqual(c.get_guess(), [2., 3., 4, 1., 0, 0.1])

        c.set_guess([6, 5, 8])
        self.assertEqual(c.get_guess(), [2., 3., 4, 6, 5, 8])

    def test_set_bounds(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        g = Gaussian()
        g2 = Gaussian()
        c.add_function(g)
        c.add_function(g2)
        lower, upper = c.get_bounds()
        self.assertEqual(lower, [0., -1., 0.0, 0, -1, 0])
        self.assertEqual(upper, [np.inf, 1., np.inf, np.inf, 1, np.inf])

        c.set_bounds([-1, -2, -3], [1, 2, 3], 0)
        lower, upper = c.get_bounds()
        self.assertEqual(lower, [-1., -2., -3, 0, -1, 0])
        self.assertEqual(upper, [1., 2, 3, np.inf, 1, np.inf])

        c.set_bounds([-4, -5, -6], [3, 4, 5])
        lower, upper = c.get_bounds()
        self.assertEqual(lower, [-1., -2., -3, -4, -5, -6])
        self.assertEqual(upper, [1., 2, 3, 3, 4, 5])


if __name__ == '__main__':
    unittest.main()
