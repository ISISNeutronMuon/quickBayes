"""Characterization tests for Four module"""
import unittest
import numpy as np
import os
from quasielasticbayes.Four import four
from quasielasticbayes.four_python import FOUR2, compress, flatten
from quasielasticbayes.fortran_python import Vec


DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class FourPythonTest(unittest.TestCase):
    """
    Characterization tests for the four python module
    just compare to the fortran results
    """

    def test_FFT_big(self):
        # reference inputs
        y = np.loadtxt(os.path.join(DATA_DIR, "FFT_test.tx"))
        y = compress(y)
        y = np.pad(y, [0, 4098-len(y)], mode="constant")
        yy = Vec(4098, True)
        yy.copy(y)
        N = 1024

        # need to convert to complex 128 as Fortran return complex 64
        tmp = four(y, N, 1, -1, -1)
        out = flatten(np.array([np.complex128(x) for x in tmp]))
        out2 = flatten(FOUR2(yy, N, 1, -1, -1))

        msd = 0
        for k in range(4098):
            msd += np.sqrt(pow(out[k]-out2[k], 2))
        max_val = np.max(out)
        tol = 0.1/100  # 0.1%
        self.assertLessEqual(msd/4098,  tol*max_val)


if __name__ == '__main__':
    unittest.main()
