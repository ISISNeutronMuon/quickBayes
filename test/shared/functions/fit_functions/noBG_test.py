import unittest
import numpy as np
from quickBayes.functions.BG import NoBG


class NoBGTest(unittest.TestCase):

    def test_no_BG_call(self):
        x = np.linspace(0, 5, 6)

        bg = NoBG()

        y = bg(x)

        for j in range(5):
            self.assertAlmostEqual(y[j], 0.0)

    def test_no_BG_report(self):
        report = {"old": [1]}

        bg = NoBG()
        out = bg.report(report)

        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 1)

    def test_N_params(self):
        bg = NoBG()
        self.assertEqual(bg.N_params, 0)

    def test_guess(self):
        bg = NoBG()
        self.assertEqual(bg.get_guess(), [])

    def test_bounds(self):
        bg = NoBG()
        lower, upper = bg.get_bounds()
        self.assertEqual(lower, [])
        self.assertEqual(upper, [])

    def test_read(self):
        report = {"old": [1]}

        bg = NoBG()
        report = bg.report(report)
        params = bg.read_from_report(report, 0)
        self.assertEqual(params, [])


if __name__ == '__main__':
    unittest.main()
