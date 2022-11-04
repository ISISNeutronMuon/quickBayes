import unittest
import numpy as np
from quasielasticbayes.v2.lorentz import Lorentzian
from quasielasticbayes.v2.BG import LinearBG
from quasielasticbayes.v2.qldata_function import QlDataFunction


class QLDataFunctionTest(unittest.TestCase):

    def test_just_bg(self):
        x = np.linspace(0, 5, 6)
        bg = LinearBG()
        ql = QlDataFunction(bg, False, x, x, 0, 6)
        y = ql(x, 1.2, 3)
        expect = 1.2*x + 3

        self.assertEqual(ql.N_params, 2)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j])

        report = {}
        report = ql.report(report, 1, 2)
        self.assertEqual(report["f1.BG gradient"], [1.])
        self.assertEqual(report["f1.BG constant"], [2.])

    def test_bg_and_delta(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        y = ql(x, 1.2, 3, .2, .1)
        expect = [-3, 0.00184, 3.076, 6.001, 9.000]

        self.assertEqual(ql.N_params, 4)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3, 4)
        self.assertEqual(report["f1.BG gradient"], [1.])
        self.assertEqual(report["f1.BG constant"], [2.])

        self.assertEqual(report["f2.f1.Amplitude"], [3.])
        self.assertEqual(report["f2.f1.Centre"], [4])

    def test_bg_and_delta_and_1_lorentzian(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)
        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, .2, .1, 1, -.2, .6)
        expect = [0.909, 0.986, 1.775, 1.076, 1.107]

        self.assertEqual(ql.N_params, 7)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3, 4, 5, 6, 7)
        self.assertEqual(report["f1.BG gradient"], [1.])
        self.assertEqual(report["f1.BG constant"], [2.])

        self.assertEqual(report["f2.f1.Amplitude"], [3.])
        self.assertEqual(report["f2.f1.Centre"], [4])

    def test_bg_and_delta_and_2_lorentzians(self):
        x = np.linspace(-5, 5, 5)
        bg = LinearBG()
        lor = Lorentzian()
        y = lor(x, 1., -.2, .6)

        ql = QlDataFunction(bg, True, x, y, -6, 6)
        ql.add_single_lorentzian()
        ql.add_single_lorentzian()

        y = ql(x, .02, 1, .2, .1, 1, -.2, .6, .7, .2, .3)
        expect = [0.913, 1.002, 2.283, 1.091, 1.111]

        self.assertEqual(ql.N_params, 10)
        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

        report = {}
        report = ql.report(report, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        self.assertEqual(report["f1.BG gradient"], [1.])
        self.assertEqual(report["f1.BG constant"], [2.])

        self.assertEqual(report["f2.f1.Amplitude"], [3.])
        self.assertEqual(report["f2.f1.Centre"], [4])

        self.assertEqual(report["f2.f2.Amplitude"], [5])
        self.assertEqual(report["f2.f2.Peak Centre"], [6])
        self.assertEqual(report["f2.f2.Gamma"], [7])

        self.assertEqual(report["f2.f3.Amplitude"], [8])
        self.assertEqual(report["f2.f3.Peak Centre"], [9])
        self.assertEqual(report["f2.f3.Gamma"], [10])


if __name__ == '__main__':
    unittest.main()
