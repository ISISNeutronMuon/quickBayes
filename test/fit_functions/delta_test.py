import unittest
import numpy as np
from quasielasticbayes.v2.delta import Delta


class DeltaTest(unittest.TestCase):

    def test_delta_call(self):
        x = np.linspace(0.0, 5.0, 11)
        d = Delta()
        y = d(x, 2.3, 1.2)
        #  x :    0, 0.5, 1,       1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.
        expect = [0, 0.0, 2.3/0.5, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0]

        for j in range(len(expect)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_delta_first_edge(self):
        x = np.linspace(0.0, 5.0, 11)
        d = Delta()
        y = d(x, 2.3, 0.2)
        #  x :    0,       0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5.
        expect = [2.3/0.5, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0]

        for j in range(len(expect)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_delta_top_edge(self):
        x = np.linspace(0.0, 5.0, 11)
        d = Delta()
        y = d(x, 2.3, 4.7)
        # this makes sure that the last bin is occupied
        #  x :    0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5,     5.
        expect = [0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 2.3/0.5, 0]

        for j in range(len(expect)):
            self.assertAlmostEqual(y[j], expect[j], 3)

    def test_delta_report(self):
        report = {"old": [1]}

        d = Delta()
        out = d.report(report, 3.2, -1)

        self.assertEqual(out["Amplitude"], [3.2])
        self.assertEqual(out["Centre"], [-1])
        self.assertEqual(out["old"], [1])
        self.assertEqual(len(out.keys()), 3)

    def test_N_params(self):
        d = Delta()
        self.assertEqual(d.N_params, 2)


if __name__ == '__main__':
    unittest.main()