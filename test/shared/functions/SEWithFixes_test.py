import unittest
import numpy as np
from quickBayes.functions.SE import StretchExp
from quickBayes.functions.SE_fix import StretchExpWithFixes


class StretchExpWithFixesTest(unittest.TestCase):

    """
    Since this function is calling a tested fit
    function, but with some fixed parameters
    we will test the values by using the original
    function. Hence, we are assuming that the
    values of the function are sufficiently
    tested by the StretchExp tests.
    Allowing us to focus on the fixes.
    """
    def test_beta(self):
        se = StretchExpWithFixes()
        self.assertEqual(se.get_beta, 0.8)

        se.set_beta(0.9)
        self.assertEqual(se.get_beta, 0.9)

    def test_init_values(self):
        se = StretchExpWithFixes(beta=0.9, FWHM=0.12)

        self.assertEqual(se.get_beta, 0.9)
        self.assertAlmostEqual(se.get_tau, 10.970, 3)

    def test_FWHM(self):
        se = StretchExpWithFixes()
        self.assertAlmostEqual(se.get_tau, 6.582, 3)

        se.set_FWHM(0.9)
        self.assertAlmostEqual(se.get_tau, 1.463, 3)

    def test_call(self):
        x = np.linspace(-0.4, 0.4, 6)

        se = StretchExp()
        se_fix = StretchExpWithFixes()

        expect = se(x, 1.0, 0.01, se.tau(.1), .7)

        se_fix.set_beta(0.7)
        se_fix.set_FWHM(.1)
        y = se_fix(x, 1.0, 0.01)

        self.assertEqual(len(y), len(expect))
        for j in range(len(y)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_report(self):
        report = {"old": [1]}

        se = StretchExpWithFixes()
        se.set_beta(0.5)
        se.set_FWHM(0.132)
        out = se.report(report, 1, 0.1)

        self.assertEqual(out["Amplitude"], [1])
        self.assertEqual(out["Peak Centre"], [0.1])
        self.assertAlmostEqual(out["tau"][0], 9.973, 3)
        self.assertEqual(out["FWHM"], [0.132])
        self.assertEqual(out["beta"], [0.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 6)

    def test_errors_report(self):
        report = {"old": [1]}

        se = StretchExpWithFixes()
        errors = np.array([0.1, 0.02])
        params = np.array([1, .2])
        out = se.report_errors(report, errors, params)

        self.assertEqual(out["Amplitude"], [0.1])
        self.assertEqual(out["Peak Centre"], [0.02])
        self.assertEqual(out["tau"], [0.])
        self.assertEqual(out["FWHM"], [0.0])
        self.assertEqual(out["beta"], [0.0])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 6)

    def test_read(self):
        report = {"old": [1]}

        se = StretchExp()
        out = se.report(report, 1, 0.1, 10, .5)

        se_fix = StretchExpWithFixes()
        params = se_fix.read_from_report(out, 0)

        self.assertEqual(params, [1, 0.1])
        self.assertEqual(se_fix.get_tau, 10)
        self.assertEqual(se_fix.get_beta, .5)

    def test_set_FWHM(self):
        se = StretchExpWithFixes(FWHM=.2)
        self.assertAlmostEqual(se.get_tau, 6.582, 3)
        se.set_FWHM(0.376)
        self.assertAlmostEqual(se.get_tau, 3.501, 3)

    def test_N_params(self):
        se = StretchExpWithFixes()
        self.assertEqual(se.N_params, 2)

    def test_guess(self):
        se = StretchExpWithFixes()
        guess = se.get_guess()
        expect = [0.1, 0.0]
        self.assertEqual(len(guess), len(expect))
        for k in range(len(guess)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

        se.set_guess([1, 2])
        self.assertEqual(se.get_guess(), [1, 2])

    def test_bounds(self):
        se = StretchExpWithFixes()
        bounds = se.get_bounds()
        self.assertEqual(bounds[0], [0, -1.])
        self.assertEqual(bounds[1], [1., 1])

        se.set_bounds([-1, -2], [3, 4])
        bounds = se.get_bounds()
        self.assertEqual(bounds[0], [-1, -2.])
        self.assertEqual(bounds[1], [3., 4])


if __name__ == '__main__':
    unittest.main()
