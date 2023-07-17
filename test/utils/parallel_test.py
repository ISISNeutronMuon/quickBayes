import unittest
from quickBayes.workflow.QSE import qse_data_main
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
    return results


class ParallelTest(unittest.TestCase):

    def test_parallel(self):

        start = time.time()
        a = []
        for j in range(6):
            a.append(function(j))
        serial = time.time() - start

        start = time.time()
        _ = parallel(6, function)
        parallel_time = time.time() - start

        self.assertLess(parallel_time, serial)


if __name__ == '__main__':
    unittest.main()
