import unittest
from quickBayes.functions.qse_fixed import QSEFixFunction
from quickBayes.functions.BG import LinearBG
from quickBayes.workflow.grid_search.qse_grid_search import QSEGridSearch
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


class QuestTest(unittest.TestCase):
    def test_quest(self):
        # get data
        sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'qse_res.npy'),
                             allow_pickle=True)
        resolution = {'x': rx, 'y': ry}

        # set values
        start_x = -0.5
        end_x = 0.5

        N_beta = 5
        beta_start = 0.75
        beta_end = 0.78

        N_FWHM = 5
        FWHM_start = 0.055
        FWHM_end = 0.056

        # create workflow
        search = QSEGridSearch()
        new_x, ry = search.preprocess_data(sx, sy, se, start_x, end_x,
                                           resolution)
        search.set_x_axis(beta_start, beta_end, N_beta, 'beta')
        search.set_y_axis(FWHM_start, FWHM_end, N_FWHM, 'FWHM')

        # create function
        bg = LinearBG()
        func = QSEFixFunction(bg, True, new_x, ry, start_x, end_x)
        func.add_single_SE()
        func.set_delta_bounds([0, -.5], [20, .5])

        # do search
        search.set_scipy_engine(func.get_guess(), *func.get_bounds())
        X, Y = search.execute(func)
        grid = search.get_grid

        # just check max value and indices
        max_val = np.max(grid)
        indices = np.where(grid == max_val)

        self.assertEqual(max_val, 1.)
        self.assertEqual(indices[0][0], 3)
        self.assertEqual(indices[1][0], 2)

        beta_slice, FWHM_slice = search.get_slices()
        expected_beta = [0.781, 0.971, 1.0, 0.904, 0.685]
        expected_FWHM = [0.814, 0.910, 0.972, 1, 0.994]

        self.assertEqual(len(beta_slice), len(expected_beta))
        self.assertEqual(len(FWHM_slice), len(expected_FWHM))
        for j in range(len(beta_slice)):
            self.assertAlmostEqual(beta_slice[j], expected_beta[j], 3)
            self.assertAlmostEqual(FWHM_slice[j], expected_FWHM[j], 3)


if __name__ == '__main__':
    unittest.main()
