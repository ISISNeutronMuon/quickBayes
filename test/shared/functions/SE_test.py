import unittest
import numpy as np
from quickBayes.functions.SE import StretchExp


class StretchExpTest(unittest.TestCase):

    def test_call_1(self):
        x = np.linspace(-0.4, 0.4, 6)

        se = StretchExp()

        y = se(x, 1.0, 0.0, 25.0, 0.5)
        # from Mantid version 6.5
        expect = [0.192, 0.299, 1.001, 1.001, 0.299, 0.192]

        self.assertEqual(len(y), len(expect))
        for j in range(len(y)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_call_2(self):
        x = np.linspace(-0.4, 0.4, 6)

        se = StretchExp()

        y = se(x, 0.5, 0.1, 25.0, 0.5)
        # from Mantid version 6.5
        expect = [0.083, 0.109, 0.201, 2.2993, 0.264, 0.122]

        self.assertEqual(len(y), len(expect))
        for j in range(len(y)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_report(self):
        report = {"old": [1]}

        se = StretchExp()
        out = se.report(report, 1, 0.1, 10, .5)

        self.assertEqual(out["Amplitude"], [1])
        self.assertEqual(out["Peak Centre"], [0.1])
        self.assertEqual(out["tau"], [10.])
        self.assertAlmostEqual(out["FWHM"][0], 0.132, 3)
        self.assertEqual(out["beta"], [0.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 6)

    def test_errors_report(self):
        report = {"old": [1]}

        se = StretchExp()
        errors = np.array([0.1, 0.02, 0.5, 0.01])
        params = np.array([1, .1, 10, .5])
        out = se.report_errors(report, errors, params)

        self.assertEqual(out["Amplitude"], [0.1])
        self.assertEqual(out["Peak Centre"], [0.02])
        self.assertEqual(out["tau"], [0.5])
        self.assertAlmostEqual(out["FWHM"][0], 0.007, 3)
        self.assertEqual(out["beta"], [0.01])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 6)

    def test_read(self):
        report = {"old": [1]}

        se = StretchExp()
        out = se.report(report, 1, 0.1, 10, .5)
        params = se.read_from_report(out, 0)

        self.assertEqual(params, [1, 0.1, 10, 0.5])

    def test_FWHM(self):
        se = StretchExp()
        FWHM = se.FWHM(3.5)
        self.assertAlmostEqual(FWHM, 0.376, 3)
        # check round trip works
        self.assertAlmostEqual(se.tau(FWHM), 3.5)

    def test_N_params(self):
        se = StretchExp()
        self.assertEqual(se.N_params, 4)

    def test_guess(self):
        se = StretchExp()
        guess = se.get_guess()
        expect = [0.1, 0.0, 6.582, 0.7]
        self.assertEqual(len(guess), len(expect))
        for k in range(len(guess)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

    def test_guess_set_guess(self):
        se = StretchExp()
        guess = se.get_guess()
        expect = [0.1, 0.0, 6.582, 0.7]
        self.assertEqual(len(guess), len(expect))
        for k in range(len(guess)):
            self.assertAlmostEqual(guess[k], expect[k], 3)
        se.set_guess([1, 2, 3, 4])
        self.assertEqual(se.get_guess(), [1, 2, 3, 4])

        se.set_guess_FWHM([.2, .1, .4, .3])
        expect = [0.2, 0.1, 3.291, 0.3]
        self.assertEqual(len(guess), len(expect))
        guess = se.get_guess()
        for k in range(len(guess)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

    def test_bounds(self):
        se = StretchExp()
        bounds = se.get_bounds()
        self.assertEqual(bounds[0], [0, -1., 0, 0])
        self.assertEqual(bounds[1], [1., 1, 100., 1.])


if __name__ == '__main__':
    unittest.main()
