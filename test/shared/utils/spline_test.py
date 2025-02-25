import unittest
from numpy import ndarray
import numpy as np
from quickBayes.utils.spline import spline


def mock_data(x: ndarray) -> ndarray:
    return 0.3*np.cos(2.*x)


class SplineTest(unittest.TestCase):

    def test_interpolation(self):
        x = np.linspace(0., 5., 30)
        y = mock_data(x)

        new_x = np.linspace(0., 5., 100)
        new_y = spline(x, y, new_x)

        expect = mock_data(new_x)

        for j in range(len(new_x)):
            self.assertAlmostEqual(new_y[j], expect[j], 3)

    def test_extrapolate(self):
        x = np.linspace(0., 5., 30)
        y = mock_data(x)

        new_x = np.linspace(-5., 10., 200)
        new_y = spline(x, y, new_x)

        expect = np.zeros(len(new_x))

        # only want non-zero's in original data range
        start = np.searchsorted(new_x, 0.)
        end = np.searchsorted(new_x, 5.)
        expect[start:end] = mock_data(new_x)[start:end]

        for j in range(len(new_x)):
            self.assertAlmostEqual(new_y[j], expect[j], 3)


if __name__ == '__main__':
    unittest.main()
