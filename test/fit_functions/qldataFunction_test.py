import unittest
import numpy as np
from quickBayes.functions.lorentz import Lorentzian
from quickBayes.functions.BG import LinearBG
from quickBayes.functions.qldata_function import QlDataFunction


class QLDataFunctionTest(unittest.TestCase):

    def test_update_resolution(self):
        def func(x):
            return x*x - 3.*x + 1.2

        x = np.linspace(-5, 5, 6)
        y = func(x)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, y, -6, 6)

        new_x = np.linspace(-5, 5, 100)
        ql.update_x_range(new_x)
        ry = ql.conv._ry
        expect = func(new_x)
        # normalise the expected values:
        expect /= sum(expect)

        self.assertEqual(len(ry), len(new_x))
        for j in range(len(ry)):
            self.assertAlmostEqual(ry[j], expect[j], 3)

    def test_get_guess(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x, 0, 6)
        with self.assertRaises(RuntimeError):
            ql.set_guess([1, 2])

    def test_just_bg(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x, 0, 6)
        y = ql(x, 1.2, 3)
        expect = 1.2*x + 3

        self.assertEqual(ql.N_params, 2)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j])

        self.assertEqual(ql.get_guess(), [0., 0.])

        bounds = ql.get_bounds()
        self.assertEqual(bounds[0], [-1, -1])
        self.assertEqual(bounds[1], [1, 1])

        report = {}
        report = ql.report(report, 1, 2)
        self.assertEqual(report["N0:f1.BG gradient"], [1.])
        self.assertEqual(report["N0:f1.BG constant"], [2.])

        ql.set_BG_guess([1., 2.])
        self.assertEqual(ql.get_guess(), [1., 2.])

    def test_bg_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        y = ql(x, 1.2, 3, .2, .1)
        expect = [-3, 0.00184, 3.076, 6.001, 9.000]

        g = ql.get_guess()
        self.assertEqual(g, [0., 0., 1., 0.])

        bounds = ql.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1])
        self.assertEqual(bounds[1], [1, 1, np.inf, 1])

        self.assertEqual(ql.N_params, 4)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3, 4)
        self.assertEqual(len(report.keys()), 4)
        self.assertEqual(report["N0:f1.BG gradient"], [1.])
        self.assertEqual(report["N0:f1.BG constant"], [2.])

        self.assertEqual(report["N0:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N0:f2.f1.Centre"], [4])

        ql.set_BG_guess([1., 2.])
        self.assertEqual(ql.get_guess(), [1., 2., 1., 0.])

        ql.set_delta_guess([3., 4.])
        self.assertEqual(ql.get_guess(), [1., 2., 3., 4.])

    def test_bg_and_delta_and_1_lorentzian(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, .2, .1, 1, .6)
        expect = [0.909, 0.987, 1.984, 1.083, 1.109]

        self.assertEqual(ql.get_guess(), [0., 0., 1., 0., 0.01, 0.02])

        bounds = ql.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0., -1, 0., 1.e-6])
        self.assertEqual(bounds[1], [1, 1, np.inf, 1, 1., 1.])

        # shared param (peak centre)
        self.assertEqual(ql.N_params, 6)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5., 6)
        self.assertEqual(len(report.keys()), 8)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N1:f2.f1.Centre"], [4])

        self.assertEqual(report["N1:f2.f2.Amplitude"], [5])
        self.assertEqual(report["N1:f2.f2.Peak Centre"], [4])
        self.assertEqual(report["N1:f2.f2.Gamma"], [6])

        self.assertEqual(report["N1:f2.f2.EISF"], [3./8.])

    def test_bg_and_1_lorentzian(self):

        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, False, x, y, -6, 6)
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, 1, 0.1, .6)
        expect = [0.909, 0.985, 1.908, 1.082, 1.108]

        self.assertEqual(ql.get_guess(), [0., 0., 0.01, 0.0, 0.02])

        bounds = ql.get_bounds()
        self.assertEqual(bounds[0], [-1, -1, 0, -1, 1.e-6])
        self.assertEqual(bounds[1], [1, 1, 1, 1, 1.])

        # shared param (peak centre)
        self.assertEqual(ql.N_params, 5)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 4, 5., 6)
        self.assertEqual(len(report.keys()), 5)
        self.assertEqual(report["N1:f1.BG gradient"], [1.])
        self.assertEqual(report["N1:f1.BG constant"], [2.])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [4])
        self.assertEqual(report["N1:f2.f1.Peak Centre"], [5])
        self.assertEqual(report["N1:f2.f1.Gamma"], [6])

    def test_read_bg_and_delta_and_1_lorentzian(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()
        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5., 6)
        params = ql.read_from_report(report, 1, 0)
        self.assertEqual(params, [1., 2., 3., 4., 5., 6.])

    def test_read_bg_and_1_lorentzian(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, False, x, y, -6, 6)
        ql.add_single_lorentzian()
        report = {}
        report = ql.report(report, 1, 2, 4, 5., 6)
        params = ql.read_from_report(report, 1, 0)
        self.assertEqual(params, [1., 2., 4., 5., 6.])

    def test_bg_and_delta_and_2_lorentzians(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)

        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, .2, .1, 1, .6, .7, .3)
        expect = [0.916, 1.016, 2.962, 1.106, 1.115]

        # shared param (peak centre)
        self.assertEqual(ql.N_params, 8)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5., 7, 8., 10)
        self.assertEqual(len(report.keys()), 12)
        self.assertEqual(report["N2:f1.BG gradient"], [1.])
        self.assertEqual(report["N2:f1.BG constant"], [2.])

        self.assertEqual(report["N2:f2.f1.Amplitude"], [3.])
        self.assertEqual(report["N2:f2.f1.Centre"], [4])

        self.assertEqual(report["N2:f2.f2.Amplitude"], [5])
        self.assertEqual(report["N2:f2.f2.Peak Centre"], [4])
        self.assertEqual(report["N2:f2.f2.Gamma"], [7])

        self.assertEqual(report["N2:f2.f3.Amplitude"], [8])
        self.assertEqual(report["N2:f2.f3.Peak Centre"], [4])
        self.assertEqual(report["N2:f2.f3.Gamma"], [10])

        self.assertEqual(report["N2:f2.f2.EISF"], [3./8.])
        self.assertEqual(report["N2:f2.f3.EISF"], [3./11.])

    def test_bg_and_2_lorentzians(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)

        ql = QlDataFunction(bg, False, x, y, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, .2, .1, 1, .6, .7, .3)
        expect = [0.907, 0.978, 1.596, 1.076, 1.107]

        # shared param (peak centre)
        self.assertEqual(ql.N_params, 7)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3., 4, 5., 6, 7)
        self.assertEqual(len(report.keys()), 8)
        self.assertEqual(report["N2:f1.BG gradient"], [1.])
        self.assertEqual(report["N2:f1.BG constant"], [2.])

        self.assertEqual(report["N2:f2.f1.Amplitude"], [3])
        self.assertEqual(report["N2:f2.f1.Peak Centre"], [4])
        self.assertEqual(report["N2:f2.f1.Gamma"], [5])

        self.assertEqual(report["N2:f2.f2.Amplitude"], [6])
        self.assertEqual(report["N2:f2.f2.Peak Centre"], [4])
        self.assertEqual(report["N2:f2.f2.Gamma"], [7])

    def test_complex_read(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        report = {}
        y = lor(x, 1., -.2, .6)

        ql = QlDataFunction(bg, True, x, y, -6, 6)
        report = ql.report(report, -1, -2, -3., -4)

        ql.add_single_lorentzian()
        report = ql.report(report, 1, 2, 3., 4, 5., 6)

        ql.add_single_lorentzian()
        report = ql.report(report, 11, 12, 13., 14, 15., 16, 17, 18)

        ql.add_single_lorentzian()
        report = ql.report(report, 21, 22, 23., 24, 25., 26, 27, 28, 29, 30)
        self.assertEqual(ql._N_peaks, 3)

        # get params for single lorentzian
        params = ql.read_from_report(report, 1, 0)
        self.assertEqual(params, [1., 2., 3., 4., 5., 6.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for three lorentzians
        params = ql.read_from_report(report, 3, 0)
        self.assertEqual(params, [21., 22., 23., 24., 25., 26., 27.,
                                  28., 29., 30.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for no lorentzians
        params = ql.read_from_report(report, 0, 0)
        self.assertEqual(params, [-1., -2., -3., -4.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for two lorentzians
        params = ql.read_from_report(report, 2, 0)
        self.assertEqual(params, [11., 12., 13., 14., 15., 16., 17., 18.])
        self.assertEqual(ql._N_peaks, 3)

    def test_complex_read_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        report = {}
        y = lor(x, 1., -.2, .6)

        ql = QlDataFunction(bg, False, x, y, -6, 6)
        report = ql.report(report, -1, -2)

        ql.add_single_lorentzian()
        report = ql.report(report, 1, 2, 4, 5., 6)

        ql.add_single_lorentzian()
        report = ql.report(report, 11, 12, 14, 15., 16, 17, 18)

        ql.add_single_lorentzian()
        report = ql.report(report, 21, 22, 24, 25., 26, 27, 28, 29, 30)
        self.assertEqual(ql._N_peaks, 3)

        # get params for single lorentzian
        params = ql.read_from_report(report, 1, 0)
        self.assertEqual(params, [1., 2., 4., 5., 6.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for three lorentzians
        params = ql.read_from_report(report, 3, 0)
        self.assertEqual(params, [21., 22., 24., 25., 26., 27., 28., 29., 30.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for no lorentzians
        params = ql.read_from_report(report, 0, 0)
        self.assertEqual(params, [-1., -2.])
        self.assertEqual(ql._N_peaks, 3)

        # get params for two lorentzians
        params = ql.read_from_report(report, 2, 0)
        self.assertEqual(params, [11., 12., 14., 15., 16., 17., 18.])
        self.assertEqual(ql._N_peaks, 3)

    def test_errors(self):
        """
        Since only the EISF is interesting for errors
        we will only test delta + 1 lorentzian
        """
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()
        report = {}
        params = [1, 2, 3, 4, 5, 6]
        errors = [.1, .2, .3, .4, .5, .6]

        report = ql.report_errors(report, errors, params)
        # EISF is the only interesting one
        self.assertEqual(len(report["N1:f2.f2.EISF"]), 1)
        self.assertAlmostEqual(report["N1:f2.f2.EISF"][0], 0.04, 2)
        # check the others
        self.assertEqual(len(report.keys()), 8)
        self.assertEqual(report["N1:f1.BG gradient"], [.1])
        self.assertEqual(report["N1:f1.BG constant"], [.2])

        self.assertEqual(report["N1:f2.f1.Amplitude"], [.3])
        self.assertEqual(report["N1:f2.f1.Centre"], [.4])

        self.assertEqual(report["N1:f2.f2.Amplitude"], [.5])
        self.assertEqual(report["N1:f2.f2.Peak Centre"], [.4])
        self.assertEqual(report["N1:f2.f2.Gamma"], [.6])

    def test_set_delta_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 3., 1])

    def test_set_delta_guess_fail(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0])

        ql.set_delta_guess([3, 1])
        self.assertEqual(ql.get_guess(), [0, 0])

    def test_set_BG_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_BG_guess([3, 1])
        self.assertEqual(ql.get_guess(), [3, 1, 1., 0])

    def test_set_func_no_peak(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 1., 0])

    def test_set_func_one_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        self.assertEqual(ql.get_guess(), [0, 0, 0.01, 0, 0.02])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 3., 2, 1])

    def test_set_func_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        self.assertEqual(ql.get_guess(), [0, 0, 1, 0, 0.01, 0.02])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 1, 2, 3., 1])

    def test_set_func_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()
        self.assertEqual(ql.get_guess(), [0, 0, 0.01, 0, 0.02, 0.01, 0.02])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 0.01, 2, 0.02, 3, 1])

        ql.set_func_guess([4, 5, -1], 0)
        self.assertEqual(ql.get_guess(), [0, 0, 4, 5, -1, 3, 1])

    def test_set_func_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()
        self.assertEqual(ql.get_guess(), [0, 0, 1, 0, 0.01, 0.02, 0.01, 0.02])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 1, 2, 0.01, 0.02, 3, 1])

        ql.set_func_guess([4, 5, -1], 0)
        self.assertEqual(ql.get_guess(), [0, 0, 1, 5, 4, -1, 3, 1])

    def test_get_func_guess(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()
        self.assertEqual(ql.get_guess(), [0, 0, 1, 0, 0.01, 0.02, 0.01, 0.02])

        ql.set_func_guess([3, 2, 1])
        self.assertEqual(ql.get_guess(), [0, 0, 1, 2, 0.01, 0.02, 3, 1])

        ql.set_func_guess([4, 5, -1], 0)
        self.assertEqual(ql.get_guess(), [0, 0, 1, 5, 4, -1, 3, 1])

        self.assertEqual(ql.get_func_guess(), [3, 5, 1])
        self.assertEqual(ql.get_func_guess(0), [4, 5, -1])

    def test_set_delta_bounds(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
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
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
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
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)

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
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
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
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_lorentzian()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 1e-6])
        self.assertEqual(upper, [1, 1, 1, 1, 1])

        ql.set_func_bounds([-3, -2, 1], [3, 2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, -3, -2, 1])
        self.assertEqual(upper, [1, 1, 3, 2, 4])

    def test_set_func_bounds_one_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_lorentzian()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 1e-6])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1, 1])

        ql.set_func_bounds([-3, -2, 1], [3, 2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -2, -3, 1])
        self.assertEqual(upper, [1, 1, np.inf, 2, 3, 4])

    def test_set_func_bounds_two_peak_no_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 1e-6, 0, 1e-6])
        self.assertEqual(upper, [1, 1, 1, 1, 1, 1, 1])

        ql.set_func_bounds([-3, -2, 1], [3, 2, 4])
        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -2, 1e-6, -3, 1])
        self.assertEqual(upper, [1, 1, 1, 2, 1, 3, 4])

    def test_set_func_bounds_two_peak_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        ql = QlDataFunction(bg, True, x, x + 1, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()

        lower, upper = ql.get_bounds()
        self.assertEqual(lower, [-1, -1, 0, -1, 0, 1e-6, 0, 1e-6])
        self.assertEqual(upper, [1, 1, np.inf, 1, 1, 1, 1, 1])

        ql.set_func_bounds([-3, -2, 1], [3, 2, 4])
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -2, 0, 1e-6, -3, 1])
        self.assertEqual(upper, [1, 1, np.inf, 2, 1, 1, 3, 4])

        ql.set_func_bounds([-5, -6, -7], [5, 6, 7], 0)
        lower, upper = ql.get_bounds()

        self.assertEqual(lower, [-1, -1, 0, -6, -5, -7, -3, 1])
        self.assertEqual(upper, [1, 1, np.inf, 6, 5, 7, 3, 4])


if __name__ == '__main__':
    unittest.main()
