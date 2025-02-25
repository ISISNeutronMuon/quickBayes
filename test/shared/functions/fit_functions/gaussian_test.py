import unittest
import numpy as np
from quickBayes.functions.gaussian import Gaussian


class GaussianTest(unittest.TestCase):

    def test_gaussian_call(self):
        x = np.linspace(0., 10., 11)

        g = Gaussian()

        y = g(x, 2.5, 4.2, 1.3)
        expect = [4e-3, 3.71e-2, 0.183, 0.501, 0.758,
                  0.635, 0.294, 7.54e-2, 1.07e-2, 8e-4, 0]

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_lorentzian_report(self):
        report = {"old": [1]}

        g = Gaussian()
        out = g.report(report, 3.2, -1, 2.5)

        self.assertEqual(out["Amplitude"], [3.2])
        self.assertEqual(out["Mean"], [-1])
        self.assertEqual(out["Sigma"], [2.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 4)

    def test_read(self):
        report = {"old": [1]}

        g = Gaussian()
        out = g.report(report, 3.2, -1, 2.5)
        params = g.read_from_report(out, 0)

        self.assertEqual(params, [3.2, -1, 2.5])

    def test_N_params(self):
        g = Gaussian()
        self.assertEqual(g.N_params, 3)

    def test_guess(self):
        g = Gaussian()
        self.assertEqual(g.get_guess(), [1.0, 0.0, 0.1])

    def test_bounds(self):
        g = Gaussian()
        bounds = g.get_bounds()
        self.assertEqual(bounds[0], [0., -1, 0.])
        self.assertEqual(bounds[1], [np.inf, 1, np.inf])

    def test_set_guess(self):
        g = Gaussian()
        self.assertEqual(g.get_guess(), [1.0, 0.0, 0.1])
        g.set_guess([2.0, 3.0, 0.4])
        self.assertEqual(g.get_guess(), [2.0, 3.0, 0.4])

    def test_set_bounds(self):
        g = Gaussian()
        bounds = g.get_bounds()
        self.assertEqual(bounds[0], [0., -1, 0.])
        self.assertEqual(bounds[1], [np.inf, 1, np.inf])

        g.set_bounds([1, 2, 3], [6, 7, 8])
        bounds = g.get_bounds()
        self.assertEqual(bounds[0], [1., 2, 3.])
        self.assertEqual(bounds[1], [6, 7, 8])


if __name__ == '__main__':
    unittest.main()
