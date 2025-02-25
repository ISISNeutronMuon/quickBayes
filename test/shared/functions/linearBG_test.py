import unittest
import numpy as np
from quickBayes.functions.BG import LinearBG


class LinearBGTest(unittest.TestCase):

    def test_linear_BG_call(self):
        x = np.linspace(0, 5, 6)

        lbg = LinearBG()

        y = lbg(x, 1.2, -3.)
        expect = [-3., -1.8, -0.6, 0.6, 1.8]

        for j in range(5):
            self.assertAlmostEqual(y[j], expect[j])

    def test_linear_BG_report(self):
        report = {"old": [1]}

        lbg = LinearBG()
        out = lbg.report(report, 3.2, -1)

        self.assertEqual(out["BG gradient"], [3.2])
        self.assertEqual(out["BG constant"], [-1])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 3)

    def test_N_params(self):
        lbg = LinearBG()
        self.assertEqual(lbg.N_params, 2)

    def test_guess(self):
        lbg = LinearBG()
        self.assertEqual(lbg.get_guess(), [0., 0.])

    def test_bounds(self):
        lbg = LinearBG()
        lower, upper = lbg.get_bounds()
        self.assertEqual(lower, [-1., -1.])
        self.assertEqual(upper, [1., 1.])

    def test_read(self):
        report = {"old": [1]}

        lbg = LinearBG()
        report = lbg.report(report, 3.2, -1)
        params = lbg.read_from_report(report, 0)
        self.assertEqual(params, [3.2, -1])

    def test_multiple_report(self):
        report = {}
        lbg = LinearBG()
        report = lbg.report(report, 3.2, -1)
        report = lbg.report(report, 4.1, .5)
        report = lbg.report(report, 0, -.4)

        self.assertEqual(report[lbg.constant], [-1, .5, -.4])
        self.assertEqual(report[lbg.grad], [3.2, 4.1, 0.])

    def test_read_index_1(self):
        report = {}
        lbg = LinearBG()
        report = lbg.report(report, 3.2, -1)
        report = lbg.report(report, 4.1, .5)
        report = lbg.report(report, 0, -.4)

        params = lbg.read_from_report(report, 1)
        self.assertEqual(params, [4.1, .5])

    def test_set_guess(self):
        lbg = LinearBG()
        self.assertEqual(lbg.get_guess(), [0., 0.])

        lbg.set_guess([1., 1.])
        self.assertEqual(lbg.get_guess(), [1., 1.])

    def test_set_bounds(self):
        lbg = LinearBG()
        lower, upper = lbg.get_bounds()
        self.assertEqual(lower, [-1., -1])
        self.assertEqual(upper, [1., 1.])

        lbg.set_bounds([0, 0], [2, 2])
        lower, upper = lbg.get_bounds()
        self.assertEqual(lower, [0., 0])
        self.assertEqual(upper, [2., 2])


if __name__ == '__main__':
    unittest.main()
