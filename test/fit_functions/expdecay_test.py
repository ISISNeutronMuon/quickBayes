import unittest
import numpy as np
from quickBayes.functions.exp_decay import ExpDecay


class ExpDecayTest(unittest.TestCase):

    def test_exp_decay_call(self):
        x = np.linspace(0., 10., 11)

        fun = ExpDecay()

        y = fun(x, 0.4, 0.2)
        expect = [0.4, 0.327, 0.268, 0.220, 0.180, 0.147,
                  0.120, 0.0986, 0.0808, 0.0661, 0.0541]

        for j in range(len(x)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_exp_decay_report(self):
        report = {"old": [1]}

        fun = ExpDecay()
        out = fun.report(report, 3.2, 2.5)

        self.assertEqual(out["Amplitude"], [3.2])
        self.assertEqual(out["lambda"], [2.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 3)

    def test_read(self):
        report = {"old": [1]}

        fun = ExpDecay()
        out = fun.report(report, 3.2, -1)
        params = fun.read_from_report(out, 0)

        self.assertEqual(params, [3.2, -1])

    def test_N_params(self):
        fun = ExpDecay()
        self.assertEqual(fun.N_params, 2)

    def test_guess(self):
        fun = ExpDecay()
        self.assertEqual(fun.get_guess(), [1.0, 0.1])

    def test_bounds(self):
        fun = ExpDecay()
        bounds = fun.get_bounds()
        self.assertEqual(bounds[0], [0., 0.001])
        self.assertEqual(bounds[1], [1., 20.])

    def test_set_guess(self):
        fun = ExpDecay()
        self.assertEqual(fun.get_guess(), [1.0, 0.1])
        fun.set_guess([2.0, 0.5])
        self.assertEqual(fun.get_guess(), [2.0, 0.5])

    def test_set_bounds(self):
        fun = ExpDecay()
        bounds = fun.get_bounds()
        self.assertEqual(bounds[0], [0., 0.001])
        self.assertEqual(bounds[1], [1., 20.])

        fun.set_bounds([1, 2], [3, 4])
        bounds = fun.get_bounds()
        self.assertEqual(bounds[0], [1, 2])
        self.assertEqual(bounds[1], [3, 4])


if __name__ == '__main__':
    unittest.main()
