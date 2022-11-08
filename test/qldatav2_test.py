import unittest
from quasielasticbayes.v2.QlData import ql_data_main
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


"""
Paramater result are from Mantid v6.5 on Windows
"""


class QlDataV2Test(unittest.TestCase):
    def test_one(self):
        sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'resolution_data_red.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        resolution = {'x': rx, 'y': ry}

        results, new_x = ql_data_main(sample, resolution,
                                      "linear", -0.4, 0.4, True)

        # not from Mantid
        self.assertAlmostEqual(results['N1:loglikelihood'][0], -659.23, 2)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -348.82, 2)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -350.36, 2)

        # from Mantid, if not then as a comment
        self.assertAlmostEqual(results['N1:f2.f2.Gamma'][0], 0.0566, 2)

        self.assertAlmostEqual(results['N2:f2.f2.Gamma'][0], 0.201, 1)  # 0.326
        self.assertAlmostEqual(results['N2:f2.f3.Gamma'][0], 0.04, 2)   # 0.046

        self.assertAlmostEqual(results['N3:f2.f2.Gamma'][0], 0.216, 2)  # 0.347
        self.assertAlmostEqual(results['N3:f2.f3.Gamma'][0], 0.050, 1)
        self.assertAlmostEqual(results['N3:f2.f4.Gamma'][0], 0.016, 2)  # 0.012

        """
        dont check the amp's directly, instead do it via EISF
        since its the ratio thats important and different defs
        have different scale factors
        """
        self.assertAlmostEqual(results['N1:f2.f2.EISF'][0], 0.0728, 2)

        self.assertAlmostEqual(results['N2:f2.f2.EISF'][0], 0.176, 1)
        self.assertAlmostEqual(results['N2:f2.f3.EISF'][0], 0.037, 3)  # 0.042

        self.assertAlmostEqual(results['N3:f2.f2.EISF'][0], 0.103, 2)  # 0.0401
        self.assertAlmostEqual(results['N3:f2.f3.EISF'][0], 0.021, 2)  # 0.0085
        self.assertAlmostEqual(results['N3:f2.f4.EISF'][0], 0.183, 2)  # 0.0668


if __name__ == '__main__':
    unittest.main()
