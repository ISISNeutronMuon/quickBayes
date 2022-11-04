import unittest
from numpy import ndarray
import numpy as np
from quasielasticbayes.v2.functions.gaussian import Gaussian
from quasielasticbayes.v2.functions.convolution import (
        ConvolutionWithResolution as conv)
from quasielasticbayes.v2.utils.crop_data import crop


def analytic(x: ndarray, amp: float, mu: float, sig: float,
             r_mu: float, r_sig: float) -> ndarray:
    """
    Expected results from wofram site on convolution
    """
    expect = np.exp(- pow(x-(mu + r_mu), 2)/(2.*(r_sig**2 + sig**2)))
    expect *= amp/np.sqrt(2.*np.pi*(r_sig**2 + sig**2))
    return expect


class ConvolutionTest(unittest.TestCase):

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

    def test_N_params(self):
        x = np.linspace(0, 1)
        c = conv(x, x, 1, 2)
        self.assertEqual(c.N_params, 0)
        g = Gaussian()
        c.add_function(g)
        self.assertEqual(c.N_params, 3)


if __name__ == '__main__':
    unittest.main()
