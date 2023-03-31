import unittest
import numpy as np
from quickBayes.functions.BG import LinearBG
from quickBayes.functions.SE import StretchExp
from quickBayes.functions.qse_fixed import QSEFixFunction


class QSEFixedFunctionTest(unittest.TestCase):

    def test_just_bg(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFixFunction(bg, False, x, x, 0, 6)
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
        qse = QSEFixFunction(bg, False, x, x, 0, 6)
        report = {}
        report = qse.report(report, 1., 2.3)
        params = qse.read_from_report(report, 0)
        self.assertEqual(params, [1., 2.3])

    def test_bg_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFixFunction(bg, True, x, y, -6, 6)
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
        qse = QSEFixFunction(bg, True, x, x, 0, 6)
        report = {}
        report = qse.report(report, 1., 2.3, .5, -.1)
        params = qse.read_from_report(report, 0)
        self.assertEqual(params, [1., 2.3, .5, -.1])

    def test_bg_and_delta_and_1_SE(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFixFunction(bg, True, x, y, -6, 6)
        qse.add_single_SE()

        y = qse(x, .02, 1, .2, .1, 1)
        expect = [0.904, 0.961, 2.445, 1.062, 1.104]

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1, 0.])
        self.assertEqual(bounds[1], [1, 1, np.inf, 1, 1])

        # shared param (peak centre)
        self.assertEqual(qse.N_params, 5)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        guess = qse.get_guess()
        expect = [0., 0., 1., 0., 0.1]
        self.assertEqual(len(guess), len(expect))
        for k in range(len(expect)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

        report = {}
        report = qse.report(report, 1, 2, 3., 4, 5.)
        self.assertEqual(len(report.keys()), 9)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N1:f2.f1.Centre"], [4])

        self.assertEqual(report["N1:f2.f2.Amplitude"], [5])
        self.assertEqual(report["N1:f2.f2.Peak Centre"], [4])
        self.assertAlmostEqual(report["N1:f2.f2.tau"][0], 6.582, 3)
        self.assertEqual(report["N1:f2.f2.FWHM"], [0.2])
        self.assertEqual(report["N1:f2.f2.beta"], [0.8])

    def test_read_bg_and_delta_and_1se(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFixFunction(bg, True, x, x, 0, 6)
        qse.add_single_SE()
        report = {}
        report = qse.report(report, 1., 2.3, .5, -.1, .1)
        params = qse.read_from_report(report, 1)
        self.assertEqual(params, [1., 2.3, .5, -.1, .1])

    def test_bg_and_1_SE(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()

        se = StretchExp()
        y = se(x, 1., 0.01, 11, .7)

        qse = QSEFixFunction(bg, False, x, y, -6, 6)
        qse.add_single_SE()

        y = qse(x, .02, 1, 1, .1)
        expect = [0.904, 0.961, 2.365, 1.062, 1.104]

        bounds = qse.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1])
        self.assertEqual(bounds[1], [1, 1, 1, 1])

        # shared param (peak centre)
        self.assertEqual(qse.N_params, 4)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        expect = [0., 0., 0.1, 0.]
        guess = qse.get_guess()
        self.assertEqual(len(guess), len(expect))
        for k in range(len(expect)):
            self.assertAlmostEqual(guess[k], expect[k], 3)

        report = {}
        report = qse.report(report, 1, 2, 3., 4)
        self.assertEqual(len(report.keys()), 7)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [3])
        self.assertEqual(report["N1:f2.f1.Peak Centre"], [4])
        self.assertAlmostEqual(report["N1:f2.f1.tau"][0], 6.582, 3)
        self.assertEqual(report["N1:f2.f1.FWHM"][0], 0.2)
        self.assertEqual(report["N1:f2.f1.beta"], [0.8])

    def test_read_bg_and_1se(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        qse = QSEFixFunction(bg, False, x, x, 0, 6)
        qse.add_single_SE()
        report = {}
        report = qse.report(report, 1., 2.3, .1, -.1)
        params = qse.read_from_report(report, 1)
        self.assertEqual(params, [1., 2.3, .1, -.1])

    def assertList(self, values, expected):
        self.assertEqual(len(values), len(expected))
        for k in range(len(values)):
            self.assertAlmostEqual(values[k], values[k], 3)

    def test_set_delta_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 3., 1])

    def test_set_delta_guess_fail(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0])

    def test_set_BG_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_BG_guess([3, 1])
        self.assertEqual(ql.get_guess(), [3, 1, 1., 0])

    def test_set_func_no_peak(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_func_guess([3, 2, 1, 4])
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

    def test_set_func_one_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 0.1, 0])

        ql.set_func_guess([3, 2])
        self.assertList(ql.get_guess(), [0, 0, 3, 2])

    def test_set_func_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1, 0, 0.1])

        ql.set_func_guess([3, 2])
        self.assertList(ql.get_guess(), [0, 0, 1, 2, 3])

    def test_set_func_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 0.1, 0, 0.1])

        ql.set_func_guess([3, 2])
        self.assertList(ql.get_guess(), [0, 0, 0.1, 2, 3])

        ql.set_func_guess([6, 5], 0)
        self.assertList(ql.get_guess(), [0, 0, 6, 5, 3])

    def test_set_func_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1., 0, 0.1, 0.1])

        ql.set_func_guess([3, 2])
        self.assertList(ql.get_guess(), [0, 0, 1., 2, 0.1, 3])

        ql.set_func_guess([4, 5], 0)
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, 3])

    def test_get_func_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()
        self.assertList(ql.get_guess(), [0, 0, 1., 0, 0.1, 0.1])

        ql.set_func_guess([3, 2])
        self.assertList(ql.get_guess(), [0, 0, 1., 2, 0.1, 3])

        ql.set_func_guess([4, 5], 0)
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, 3])
        self.assertList(ql.get_guess(), [0, 0, 1., 5, 4, 3])

        self.assertEqual(ql.get_func_guess(), [3, 5])
        self.assertEqual(ql.get_func_guess(0), [4, 5])

    def test_set_delta_bounds(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
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
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
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
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)

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
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
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
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1])
        self.assertEqual(upper, [1, 1, 1, 1])

        ql.set_func_bounds([-3, -2], [3, 2])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, -3, -2])
        self.assertEqual(upper, [1, 1, 3, 2])

    def test_set_func_bounds_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1])

        ql.set_func_bounds([-3, -2], [3, 2])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -2, -3])
        self.assertEqual(upper, [1, 1, np.inf, 2, 3])

    def test_set_func_bounds_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0])
        self.assertEqual(upper, [1, 1, 1, 1, 1])

        ql.set_func_bounds([-3, -2], [3, 2])
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -2, -3])
        self.assertEqual(upper, [1, 1, 1, 2, 3])

    def test_set_func_bounds_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 0])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1, 1])

        ql.set_func_bounds([-3, -2], [3, 2])
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -2, 0, -3])
        self.assertEqual(upper, [1, 1, np.inf, 2, 1, 3])

        ql.set_func_bounds([-5, -6], [5, 6], 0)
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -6, -5, -3])
        self.assertEqual(upper, [1, 1, np.inf, 6, 5, 3])

    @staticmethod
    def get_se(func, index):
        return func.conv._funcs[index]

    def test_beta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()

        self.assertEqual(self.get_se(ql, 1)._beta, 0.8)

        ql.set_beta(0.9)

        self.assertEqual(self.get_se(ql, 1)._beta, 0.9)

        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5.)
        self.assertEqual(report["N1:f2.f2.beta"], [0.9])

        ql.set_beta(1.0)
        self.assertEqual(self.get_se(ql, 1)._beta, 1.0)

        _ = ql.read_from_report(report, 1)

        self.assertEqual(self.get_se(ql, 1)._beta, 0.9)

    def test_2_betas(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        self.assertEqual(self.get_se(ql, 1).get_beta, 0.8)
        self.assertEqual(self.get_se(ql, 2).get_beta, 0.8)

        ql.set_beta(0.9)
        self.assertEqual(self.get_se(ql, 1).get_beta, 0.8)
        self.assertEqual(self.get_se(ql, 2).get_beta, 0.9)

        ql.set_beta(1.0, 0)
        self.assertEqual(self.get_se(ql, 1).get_beta, 1.0)
        self.assertEqual(self.get_se(ql, 2).get_beta, 0.9)

    def test_tau(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()

        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 6.582, 3)

        ql.set_FWHM(0.4)

        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 3.291, 3)

        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5.)
        self.assertEqual(report["N1:f2.f2.FWHM"], [0.4])

        ql.set_FWHM(.1)
        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 13.164, 3)

        _ = ql.read_from_report(report, 1)

        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 3.291, 3)

    def test_2_taus(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QSEFixFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_SE()
        ql.add_single_SE()

        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 6.582, 3)
        self.assertAlmostEqual(self.get_se(ql, 2).get_tau, 6.582, 3)

        ql.set_FWHM(0.4)

        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 6.582, 3)
        self.assertAlmostEqual(self.get_se(ql, 2).get_tau, 3.291, 3)

        ql.set_FWHM(.1, 0)
        self.assertAlmostEqual(self.get_se(ql, 1).get_tau, 13.164, 3)
        self.assertAlmostEqual(self.get_se(ql, 2).get_tau, 3.291, 3)


if __name__ == '__main__':
    unittest.main()
