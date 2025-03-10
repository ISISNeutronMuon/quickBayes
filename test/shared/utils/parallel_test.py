import unittest
from quickBayes.workflow.model_selection.QSE import qse_data_main
from quickBayes.utils.parallel import parallel
import numpy as np
import os.path
import time


DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


def function(j):

    sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
    rx, ry, re = np.load(os.path.join(DATA_DIR, 'qse_res.npy'),
                         allow_pickle=True)

    sample = {'x': sx, 'y': sy, 'e': se}
    resolution = {'x': rx, 'y': ry}
    results = {}
    errors = {}
    (results, errors,
     new_x, fits, fit_e) = qse_data_main(sample, resolution,
                                         "linear", -0.4, 0.4, True,
                                         results, errors)
    return results, j


class ParallelTest(unittest.TestCase):

    def test_parallelSpeed(self):

        start = time.time()
        a = []
        for j in range(6):
            a.append(function(j))
        serial = time.time() - start

        start = time.time()
        _ = parallel(list(range(6)), function)
        parallel_time = time.time() - start

        self.assertLess(parallel_time, serial)

    def test_parallelResults(self):

        data = parallel(list(range(2)), function)
        first = data[0][0]
        second = data[1][0]
        self.assertEqual(first.keys(), second.keys())
        for key in first.keys():
            self.assertEqual(first[key], second[key])

        # check the j's
        self.assertEqual(data[0][1], 0)
        self.assertEqual(data[1][1], 1)


if __name__ == '__main__':
    unittest.main()
