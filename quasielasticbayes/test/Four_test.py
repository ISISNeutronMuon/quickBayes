"""Characterization tests for Four module"""
import unittest
import numpy as np
from quasielasticbayes.Four import four
from quasielasticbayes.python.four import FOUR2, compress, flatten
from quasielasticbayes.python.fortran_python import *
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class FourTest(unittest.TestCase):
    """
    Characterization tests for the four Fortran module
    """
    # P = +1, N = -1, Z = 0
    def test_Four_PPN_Re(self):
        # reference inputs
        x = np.zeros(4098, dtype=complex)
        nonzero_entries = 8
        x[:nonzero_entries] = [complex(np.cos(float(i)+0.1), 0.)
                               for i in range(nonzero_entries)]
        y = x
        out = four(y, 8, 1, 1, -1)
        dp = 3
        self.assertAlmostEqual(-1.6806+3.6243j, out[0], dp)
        self.assertAlmostEqual(1.4300-0.4846j, out[1], dp)
        self.assertAlmostEqual(0.5016-0.4846j, out[2], dp)
        self.assertAlmostEqual(1.4299+3.6243j, out[3], dp)
        self.assertAlmostEqual(-0.5748, out[4], dp)
        self.assertAlmostEqual(0.3780, out[5], dp)
        self.assertAlmostEqual(0.9833, out[6], dp)
        self.assertAlmostEqual(0.6845, out[7], dp)

    # In the code it only seems to use 1,-1,-1 and 1,1,0
    # So these are the two we will test
    def test_Four_PNN_Re(self):
        x = np.zeros(4098, dtype=complex)
        nonzero_entries = 8
        x[:nonzero_entries] = [complex(np.cos(float(i)+0.1), 0.)
                               for i in range(nonzero_entries)]
        y = x
        out = four(y, 8, 1, -1, -1)
        dp = 3
        self.assertAlmostEqual(-1.681+3.6243j, out[0], dp)
        self.assertAlmostEqual(1.4299-0.4846j, out[1], dp)
        self.assertAlmostEqual(0.5016-0.4846j, out[2], dp)
        self.assertAlmostEqual(1.4299+3.6243j, out[3], dp)
        self.assertAlmostEqual(-0.5748, out[4], dp)
        self.assertAlmostEqual(0.37800, out[5], dp)
        self.assertAlmostEqual(0.9833, out[6], dp)
        self.assertAlmostEqual(0.6845, out[7], dp)

    def test_Four_PNN_Im(self):
        x = np.zeros(4098, dtype=complex)
        nonzero_entries = 8
        x[:nonzero_entries] = [complex(np.cos(float(i)+0.1), i)
                               for i in range(nonzero_entries)]
        y = x
        out = four(y, 8, 1, -1, -1)
        dp = 3
        self.assertAlmostEqual(-1.681+13.2812j, out[0], dp)
        self.assertAlmostEqual(-2.5701+1.1722j, out[1], dp)
        self.assertAlmostEqual(0.5016-2.1415j, out[2], dp)
        self.assertAlmostEqual(5.4299-6.0326j, out[3], dp)
        self.assertAlmostEqual(-0.5748+4j, out[4], dp)
        self.assertAlmostEqual(0.3780+5j, out[5], dp)
        self.assertAlmostEqual(0.9833+6j, out[6], dp)
        self.assertAlmostEqual(0.6845+7j, out[7], dp)

    def test_Four_PPZ_Re(self):
        x = np.zeros(4098, dtype=complex)
        nonzero_entries = 8
        x[:nonzero_entries] = [complex(np.cos(float(i)+0.1), 0.)
                               for i in range(nonzero_entries)]
        y = x
        out = four(y, 8, 1, 1, 0)
        dp = 3
        self.assertAlmostEqual(-0.0553, out[0], dp)
        self.assertAlmostEqual(1.5000+1.4527j, out[1], dp)
        self.assertAlmostEqual(1.0357, out[2], dp)
        self.assertAlmostEqual(1.5000-1.4527j, out[3], dp)
        self.assertAlmostEqual(-0.0554, out[4], dp)
        self.assertAlmostEqual(0.3780, out[5], dp)
        self.assertAlmostEqual(0.9833, out[6], dp)
        self.assertAlmostEqual(0.6845, out[7], dp)

    def test_Four_PPZ_Im(self):
        x = np.zeros(4098, dtype=complex)
        nonzero_entries = 8
        x[:nonzero_entries] = [complex(np.cos(float(i)+0.1), i)
                               for i in range(nonzero_entries)]
        y = x
        out = four(y, 8, 1, 1, 0)
        dp = 3
        self.assertAlmostEqual(5.9446, out[0], dp)
        self.assertAlmostEqual(1.4999-1.3757j, out[1], dp)
        self.assertAlmostEqual(1.0357-2j, out[2], dp)
        self.assertAlmostEqual(1.4999-4.2812j, out[3], dp)
        self.assertAlmostEqual(-6.0554, out[4], dp)
        self.assertAlmostEqual(0.3780+5j, out[5], dp)
        self.assertAlmostEqual(0.9833+6j, out[6], dp)
        self.assertAlmostEqual(0.6845+7j, out[7], dp)	

    def test_FFT_big(self):
        # reference inputs
        y = np.loadtxt("C:\\Users\\BTR75544\\work\\quasielasticbayes\\quasielasticbayes\\test\\data\\FFT_test.tx")
        y = compress(y)
        y=np.pad(y,[0, 4098-len(y)], mode="constant")
        yy = vec(4098, True)
        yy.copy(y)
        N = 1024

        # need to convert to complex 128 as Fortran return complex 64
        tmp = four(y, N, 1, -1, -1)
        out = flatten(np.array([np.complex128(x) for x in tmp]))
        out2 = flatten(FOUR2(yy,N,1,-1,-1))

        x = np.asarray([k for k in range(len(out))])
        x2 = np.asarray([k for k in range(len(out2))])

        s=0
        msd = 0
        for k in range(4098):
            msd += np.sqrt(pow(out[k]-out2[k],2))
        max_val = np.max(out)
        tol = 0.1/100 # 0.1%
        self.assertLessEqual(msd/4098,  tol*max_val)

if __name__ == '__main__':
    unittest.main()
