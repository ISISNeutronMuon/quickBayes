"""Characterization tests for QLres module"""
import os.path
import unittest
import numpy as np
from quasielasticbayes.testing import load_json

from quasielasticbayes.Four import four
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class FFTTest(unittest.TestCase):
    """
    Characterization tests using inputs that have been accepted as correct.
    The output is based on running BayesQuasi algorithm in mantid 6.2
    with the inputs taken from the BayesQuasiTest unit test
    """

    def test_FFT(self):
        # reference inputs
        x = np.asarray([k for k in range(4)])
        y = np.cos(x)
        out = four(y,4,1,1,1)
        print(out)

if __name__ == '__main__':
    unittest.main()