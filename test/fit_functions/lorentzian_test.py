import unittest
import numpy as np
from quasielasticbayes.v2.lorentz import Lorentzian


class LorentzianTest(unittest.TestCase):

    def test_lorentzian_call(self):
        x = np.linspace(-0.4, 0.4, 6)

        lor = Lorentzian()

        y = lor(x, 20.3, 0.031, 0.3)
        # from Mantid version 6.5
        expect = [4.654, 10.103, 27.835, 38.924, 14.645, 6.109]

        for j in range(5):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_lorentzian_report(self):
        report = {"old": [1]}

        lor = Lorentzian()
        out = lor.report(report, 3.2, -1, 2.5)

        self.assertEqual(out["Amplitude"], [3.2])
        self.assertEqual(out["Peak Centre"], [-1])
        self.assertEqual(out["Gamma"], [2.5])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 4)

    def test_N_params(self):
        lor = Lorentzian()
        self.assertEqual(lor.N_params, 3)


if __name__ == '__main__':
    unittest.main()
