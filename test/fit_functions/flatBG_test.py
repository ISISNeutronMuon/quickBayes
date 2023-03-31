import unittest
import numpy as np
from quickBayes.functions.BG import FlatBG


class FlatBGTest(unittest.TestCase):

    def test_flat_BG_call(self):
        x = np.linspace(0, 5, 6)

        bg = FlatBG()

        y = bg(x, 1.2)

        for j in range(5):
            self.assertAlmostEqual(y[j], 1.2)

    def test_flat_BG_report(self):
        report = {"old": [1]}

        bg = FlatBG()
        out = bg.report(report, 3.2)

        self.assertEqual(out["BG constant"], [3.2])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 2)

    def test_flat_BG_report_error(self):
        report = {"old": [1]}

        bg = FlatBG()
        out = bg.report_errors(report, [0.1], [3.2])

        self.assertEqual(out["BG constant"], [0.1])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 2)

    def test_N_params(self):
        bg = FlatBG()
        self.assertEqual(bg.N_params, 1)

    def test_guess(self):
        bg = FlatBG()
        self.assertEqual(bg.get_guess(), [0.])

    def test_bounds(self):
        bg = FlatBG()
        lower, upper = bg.get_bounds()
        self.assertEqual(lower, [-1.])
        self.assertEqual(upper, [1.])

    def test_read(self):
        report = {"old": [1]}

        bg = FlatBG()
        report = bg.report(report, 3.2)
        params = bg.read_from_report(report, 0)
        self.assertEqual(params, [3.2])

    def test_set_guess(self):
        bg = FlatBG()
        self.assertEqual(bg.get_guess(), [0.])

        bg.set_guess([1.])
        self.assertEqual(bg.get_guess(), [1.])

    def test_set_bounds(self):
        bg = FlatBG()
        lower, upper = bg.get_bounds()
        self.assertEqual(lower, [-1.])
        self.assertEqual(upper, [1.])

        bg.set_bounds([0], [2])
        lower, upper = bg.get_bounds()
        self.assertEqual(lower, [0.])
        self.assertEqual(upper, [2.])


if __name__ == '__main__':
    unittest.main()
