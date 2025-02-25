import unittest
import numpy as np
from quickBayes.functions.SE import StretchExp
from quickBayes.functions.lorentz import Lorentzian
from quickBayes.functions.BG import LinearBG
from quickBayes.functions.composite import CompositeFunction


class CompositeFunctionTest(unittest.TestCase):

    def test_empty(self):
        x = np.linspace(0, 5, 6)
        c = CompositeFunction()
        y = c(x)

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], 0.0, 8)

    def test_add_function(self):
        x = np.linspace(0, 5, 6)
        c = CompositeFunction()
        bg = LinearBG('f0.')

        self.assertEqual(c.N_params, 0)
        c.add_function(bg)
        self.assertEqual(c.N_params, 2)

        a = 1.3
        b = -3.1
        expect = bg(x, a, b)
        y = c(x, a, b)

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 8)

    def test_too_few_params(self):
        x = np.linspace(0, 5, 6)
        c = CompositeFunction()
        bg = LinearBG('f0.')
        c.add_function(bg)
        with self.assertRaises(ValueError):
            _ = c(x, 1.)  # should have 2 params

    def test_too_many_params(self):
        x = np.linspace(0, 5, 6)
        c = CompositeFunction()
        bg = LinearBG('f0.')
        c.add_function(bg)

        with self.assertRaises(ValueError):
            _ = c(x, 1, 2, 3)  # should have 2 params

    def test_sum(self):
        x = np.linspace(-0.4, 0.4, 6)

        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG('f0.')
        c.add_function(bg)
        c.add_function(lor)

        y_l = lor(x, 20.3, 0.031, 0.3)
        y_bg = bg(x, 1, -2)

        y = c(x, 1, -2, 20.3, 0.031, 0.3)

        for j in range(5):
            self.assertAlmostEqual(y[j], y_l[j] + y_bg[j], 3)

    def test_report(self):
        report = {"old": [1]}

        lor = Lorentzian()
        c = CompositeFunction()
        c.add_function(lor)

        out = c.report(report, 3.2, -1, 2.5)
        self.assertEqual(out["f1.Amplitude"], [3.2])
        self.assertEqual(out["f1.Peak Centre"], [-1])
        self.assertEqual(out["f1.Gamma"], [2.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 4)

    def test_report2(self):
        report = {"old": [1]}

        lor = Lorentzian()
        c = CompositeFunction()
        c.add_function(lor)
        lor2 = Lorentzian()
        c.add_function(lor2)

        out = c.report(report, 3.2, -1, 2.5, 5, 6, 7)

        self.assertEqual(out["f1.Amplitude"], [3.2])
        self.assertEqual(out["f1.Peak Centre"], [-1])
        self.assertEqual(out["f1.Gamma"], [2.5])

        self.assertEqual(out["f2.Amplitude"], [5])
        self.assertEqual(out["f2.Peak Centre"], [6])
        self.assertEqual(out["f2.Gamma"], [7])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 7)

    def test_custom_report_errors(self):
        """
        stretch exp calculates the error for FWHM,
        test that this is used when its part of
        a composite function.
        """
        se = StretchExp()
        c = CompositeFunction()
        c.add_function(se)

        params = [1.97e-1, -1.43e-3, 2.10e1, 7.73e-1]
        sigma = [3.e-4, 5.e-5, 5.e-2, 1.7e-3]

        errors = c.report_errors({}, sigma, params)

        self.assertAlmostEqual(errors['f1.FWHM'][0], 0.00015, 5)

    def test_read(self):
        report = {"old": [1]}

        lor = Lorentzian()
        c = CompositeFunction()
        c.add_function(lor)

        out = c.report(report, 3.2, -1, 2.5)
        params = c.read_from_report(out, 0)
        self.assertEqual(params, [3.2, -1, 2.5])

    def test_read2(self):
        report = {"old": [1]}

        lor = Lorentzian()
        c = CompositeFunction()
        c.add_function(lor)
        lor2 = Lorentzian()
        c.add_function(lor2)

        out = c.report(report, 3.2, -1, 2.5, 5, 6, 7)
        params = c.read_from_report(out, 0)
        self.assertEqual(params, [3.2, -1, 2.5, 5, 6, 7])

    def test_guess(self):
        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG()
        c.add_function(bg)
        c.add_function(lor)

        self.assertEqual(c.get_guess(), [0., 0., 0.01, 0., 0.02])

    def test_bounds(self):
        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG()
        c.add_function(bg)
        c.add_function(lor)

        bounds = c.get_bounds()

        self.assertEqual(bounds[0], [-1., -1., 0., -1., 1.e-6])
        self.assertEqual(bounds[1], [1., 1., 1., 1., 1.])

    def test_set_guess(self):
        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG()
        c.add_function(bg)
        c.add_function(lor)

        self.assertEqual(c.get_guess(), [0., 0., 0.01, 0., 0.02])

        c.set_guess([1, 2], 0)
        self.assertEqual(c.get_guess(), [1., 2., 0.01, 0., 0.02])

        c.set_guess([3, 4, 5])
        self.assertEqual(c.get_guess(), [1., 2., 3, 4., 5])

    def test_set_bounds(self):
        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG()
        c.add_function(bg)
        c.add_function(lor)

        bounds = c.get_bounds()

        self.assertEqual(bounds[0], [-1., -1., 0., -1., 1.e-6])
        self.assertEqual(bounds[1], [1., 1., 1., 1., 1.])

        c.set_bounds([-2, -3], [2, 3], 0)
        bounds = c.get_bounds()

        self.assertEqual(bounds[0], [-2., -3., 0., -1., 1.e-6])
        self.assertEqual(bounds[1], [2., 3., 1., 1., 1.])

        c.set_bounds([-4, -5, -6], [4, 5, 6])
        bounds = c.get_bounds()

        self.assertEqual(bounds[0], [-2., -3., -4., -5., -6])
        self.assertEqual(bounds[1], [2., 3., 4., 5., 6.])

    def test_update_prefix(self):
        lor = Lorentzian()
        c = CompositeFunction()
        bg = LinearBG()
        c.add_function(bg)
        c.add_function(lor)

        for j, fun in enumerate(c._funcs):
            self.assertEqual(fun._prefix, f'f{j+1}.')

        c.update_prefix("test:")

        for j, fun in enumerate(c._funcs):
            self.assertEqual(fun._prefix, f'test:f{j+1}.')


if __name__ == '__main__':
    unittest.main()
