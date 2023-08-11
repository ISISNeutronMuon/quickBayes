import unittest
import numpy as np
from quickBayes.functions.BG import LinearBG
from quickBayes.functions.SE import StretchExp
from quickBayes.functions.qse_function import QSEFunction


class QSEFunctionTest(unittest.TestCase):

    def test_just_bg(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFunction(bg, False, x, x, 0, 6)
        y = qse(x, 1.2, 3)
        expect = 1.2*x + 3

        self.assertEqual(qse.N_params, 2)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j])

        self.assertEqual(qse.get_guess(), [0., 0.])

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1])
        self.assertEqual(bounds[1], [1, 1])

        report = {}
        report = qse.report(report, 1, 2)
        self.assertEqual(report["N0:f1.BG gradient"], [1.])
        self.assertEqual(report["N0:f1.BG constant"], [2.])

    def test_read_just_bg(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFunction(bg, False, x, x, 0, 6)
        report = {}
        report = qse.report(report, 1., 2.3)
        params = qse.read_from_report(report, 0)
        self.assertEqual(params, [1., 2.3])

    def test_bg_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFunction(bg, True, x, y, -6, 6)
        y = qse(x, 1.2, 3, .2, .1)
        expect = [-3, 0.0001, 3.080, 6.000, 9.000]

        self.assertEqual(qse.get_guess(), [0., 0., 1., 0.])

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1])
        self.assertEqual(bounds[1], [1, 1, np.inf, 1])

        self.assertEqual(qse.N_params, 4)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = qse.report(report, 1, 2, 3, 4)
        self.assertEqual(len(report.keys()), 4)
        self.assertEqual(report["N0:f1.BG gradient"], [1.])
        self.assertEqual(report["N0:f1.BG constant"], [2.])

        self.assertEqual(report["N0:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N0:f2.f1.Centre"], [4])

    def test_read_bg_and_delta(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFunction(bg, True, x, x, 0, 6)
        report = {}
        report = qse.report(report, 1., 2.3, .5, -.1)
        params = qse.read_from_report(report, 0)
        self.assertEqual(params, [1., 2.3, .5, -.1])

    def test_bg_and_delta_and_1_SE(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFunction(bg, True, x, y, -6, 6)
        qse.add_single_SE()

        y = qse(x, .02, 1, .2, .1, 1, 10., 0.5)
        expect = [0.911, 0.972, 2.345, 1.074, 1.112]

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1, 0., 0, 0])
        self.assertEqual(bounds[1], [1, 1, np.inf, 1, 1, 100, 1])

        # shared param (peak centre)
        self.assertEqual(qse.N_params, 7)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        guess = qse.get_guess()
        expect = [0., 0., 1., 0., 0.1, 6.582, 0.7]
        self.assertEqual(len(guess), len(expect))
        for k in range(len(expect)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

        report = {}
        report = qse.report(report, 1, 2, 3., 4, 5., 6, 7)
        self.assertEqual(len(report.keys()), 9)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N1:f2.f1.Centre"], [4])

        self.assertEqual(report["N1:f2.f2.Amplitude"], [5])
        self.assertEqual(report["N1:f2.f2.Peak Centre"], [4])
        self.assertEqual(report["N1:f2.f2.tau"], [6])
        self.assertAlmostEqual(report["N1:f2.f2.FWHM"][0], 0.219, 3)
        self.assertEqual(report["N1:f2.f2.beta"], [7])

    def test_report_errors(self):
        """
        stretch exp calculates the error for FWHM,
        test that this is used when its part of
        a quasielastic function.
        """

        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFunction(bg, True, x, y, -6, 6)
        qse.add_single_SE()

        y = qse(x, .02, 1, .2, .1, 1, 10., 0.5)

        params = [0., 0., 0.911, 0.972, 2.345, 2.10e1, 7.73e-1]
        sigma = [0., 0., 1.2e-2, 3.e-4, 5.e-5, 5.e-2, 1.7e-3]

        errors = qse.report_errors({}, sigma, params)
        self.assertAlmostEqual(errors['N1:f2.f2.FWHM'][0], 0.00015, 5)

    def test_read_bg_and_delta_and_1se(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFunction(bg, True, x, x, 0, 6)
        qse.add_single_SE()
        report = {}
        report = qse.report(report, 1., 2.3, .5, -.1, .1, 10, .7)
        params = qse.read_from_report(report, 1)
        self.assertEqual(params, [1., 2.3, .5, -.1, .1, 10, .7])

    def test_bg_and_1_SE(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFunction(bg, False, x, y, -6, 6)
        qse.add_single_SE()

        y = qse(x, .02, 1, 1, .1, 10., 0.5)
        expect = [0.911, 0.972, 2.265, 1.074, 1.112]

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1, 0, 0])
        self.assertEqual(bounds[1], [1, 1, 1, 1, 100, 1])

        # shared param (peak centre)
        self.assertEqual(qse.N_params, 6)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        expect = [0., 0., 0.1, 0., 6.582, 0.7]
        guess = qse.get_guess()
        self.assertEqual(len(guess), len(expect))
        for k in range(len(expect)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

        report = {}
        report = qse.report(report, 1, 2, 3., 4, 5., 6)
        self.assertEqual(len(report.keys()), 7)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [3])
        self.assertEqual(report["N1:f2.f1.Peak Centre"], [4])
        self.assertEqual(report["N1:f2.f1.tau"], [5])
        self.assertAlmostEqual(report["N1:f2.f1.FWHM"][0], 0.263, 3)
        self.assertEqual(report["N1:f2.f1.beta"], [6])

    def test_read_bg_and_1se(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFunction(bg, False, x, x, 0, 6)
        qse.add_single_SE()
        report = {}
        report = qse.report(report, 1., 2.3, .1, -.1, 10, .7)
        params = qse.read_from_report(report, 1)
        self.assertEqual(params, [1., 2.3, .1, -.1, 10, .7])

    def assertList(self, values, expected):
        self.assertEqual(len(values), len(expected))
        for k in range(len(values)):
            self.assertAlmostEqual(values[k], values[k], 3)

    def test_set_delta_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 3., 1])

    def test_set_delta_guess_fail(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0])

    def test_set_BG_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_BG_guess([3, 1])
        self.assertEqual(ql.get_guess(), [3, 1, 1., 0])

    def test_set_func_no_peak(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

    def test_set_func_one_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 0.1, 0, 6.582, 0.7])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 3, 2, 1, 4])

    def test_set_func_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1, 0, 0.1, 6.582, 0.7])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 1, 2, 3, 1, 4])

    def test_set_func_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 0.1, 0, 6.582, 0.7,
                                         0.1, 6.582, 0.7])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 0.1, 2, 6.582, 0.7, 3, 1, 4])

        ql.set_func_guess([6, 5, -1, -2], 0)
        self.assertList(ql.get_guess(), [0, 0, 6, 5, -1, -2, 3, 1, 4])

    def test_set_func_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1., 0, 0.1,
                                         6.582, 0.7, 0.1, 6.582, 0.7])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 1., 2, 0.1,
                                         6.582, 0.7, 3, 1, 4])

        ql.set_func_guess([4, 5, -1, -2], 0)
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, -1, -2, 3, 1, 4])

    def test_set_func_guess_FWHM(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1, 0, 0.1, 6.582, 0.7])
        ql.set_func_guess_FWHM([3, 2, 0.4, 4])
        self.assertList(ql.get_guess(), [0, 0, 1, 2, 3, 3.291, 4])

    def test_get_func_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1., 0, 0.1,
                                         6.582, 0.7, 0.1, 6.582, 0.7])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 1., 2, 0.1,
                                         6.582, 0.7, 3, 1, 4])

        ql.set_func_guess([4, 5, -1, -2], 0)
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, -1, -2, 3, 1, 4])
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, -1, -2, 3, 1, 4])

        self.assertEqual(ql.get_func_guess(), [3, 5, 1, 4])
        self.assertEqual(ql.get_func_guess(0), [4, 5, -1, -2])

    def test_set_delta_bounds(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1])
        self.assertEqual(upper, [1, 1, np.inf, 1])

        ql.set_delta_bounds([-3, -2], [2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, -3, -2])
        self.assertEqual(upper, [1, 1, 2, 4])

    def test_set_delta_bounds_fail(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1])
        self.assertEqual(upper, [1, 1])

        ql.set_delta_bounds([-3, -2], [2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1])
        self.assertEqual(upper, [1, 1])

    def test_set_BG_bounds(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1])
        self.assertEqual(upper, [1, 1, np.inf, 1])

        ql.set_BG_bounds([-3, -2], [2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-3, -2, 0, -1])
        self.assertEqual(upper, [2, 4, np.inf, 1])

    def test_set_func_bounds_no_peak(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1])
        self.assertEqual(upper, [1, 1, np.inf, 1])

        ql.set_func_bounds([-3, -2, 1], [1, 2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1])
        self.assertEqual(upper, [1, 1, np.inf, 1])

    def test_set_func_bounds_one_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 0])
        self.assertEqual(upper, [1, 1, 1, 1, 100, 1])

        ql.set_func_bounds([-3, -2, 1, -4], [3, 2, 4, 5])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, -3, -2, 1, -4])
        self.assertEqual(upper, [1, 1, 3, 2, 4, 5])

    def test_set_func_bounds_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 0, 0])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1, 100, 1])

        ql.set_func_bounds([-3, -2, 1, -4], [3, 2, 4, 5])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -2, -3, 1, -4])
        self.assertEqual(upper, [1, 1, np.inf, 2, 3, 4, 5])

    def test_set_func_bounds_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 0, 0, 0, 0])
        self.assertEqual(upper, [1, 1, 1, 1, 100, 1, 1, 100, 1])

        ql.set_func_bounds([-3, -2, 1, -4], [3, 2, 4, 5])
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -2, 0, 0, -3, 1, -4])
        self.assertEqual(upper, [1, 1, 1, 2, 100, 1, 3, 4, 5])

    def test_set_func_bounds_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 0, 0, 0, 0, 0])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1, 100, 1, 1, 100, 1])

        ql.set_func_bounds([-3, -2, 1, -4], [3, 2, 4, 5])
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -2, 0, 0, 0, -3, 1, -4])
        self.assertEqual(upper, [1, 1, np.inf, 2, 1, 100, 1, 3, 4, 5])

        ql.set_func_bounds([-5, -6, -7, -8], [5, 6, 7, 8], 0)
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -6, -5, -7, -8, -3, 1, -4])
        self.assertEqual(upper, [1, 1, np.inf, 6, 5, 7, 8, 3, 4, 5])


if __name__ == '__main__':
    unittest.main()
