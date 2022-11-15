import unittest
import numpy as np
from quasielasticbayes.v2.functions.BG import LinearBG


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


if __name__ == '__main__':
    unittest.main()
