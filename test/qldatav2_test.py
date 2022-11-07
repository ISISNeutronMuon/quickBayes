import unittest
from quasielasticbayes.v2.QlData import ql_data_main
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


"""
All result are from Mantid v6.5 on Windows
"""


class QlDataV2Test(unittest.TestCase):
    def test_one(self):
        sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'resolution_data_red.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        resolution = {'x': rx, 'y': ry}

        results, probs = ql_data_main(sample, resolution,
                                      "linear", -0.4, 0.4, True)

        self.assertAlmostEqual(probs[0], -68736, 1)
        self.assertAlomostEqual(list(results.keys()), [])


if __name__ == '__main__':
    unittest.main()
