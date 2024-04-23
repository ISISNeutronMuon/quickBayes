import unittest
from quickBayes.functions.qse_function import QSEFunction
from quickBayes.functions.BG import LinearBG
from quickBayes.workflow.model_selection.QSE import (qse_data_main,
                                                     QlStretchedExp)
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


"""
Paramater result are from Mantid v6.5 on Windows
"""


class QSETest(unittest.TestCase):

    def test_raw(self):
        sx = [0, 1, 3, 5, 6.9]
        sy = [1, 2, .2, .4, .1]
        se = [.1, .2, .3, .4, .5]

        sample = {'x': sx, 'y': sy, 'e': se}
        res = {'x': sx, 'y': sy}

        results = {}
        errors = {}
        workflow = QlStretchedExp(results, errors)
        new_x, ry = workflow.preprocess_data(sample['x'], sample['y'],
                                             sample['e'],
                                             0., 7., res)
        expect_x = [0, 1.17, 2.33, 3.5, 4.67, 5.83, 7]
        self.assertEqual(len(new_x), len(expect_x))
        for k in range(len(new_x)):
            self.assertAlmostEqual(new_x[k], expect_x[k], 2)

        raw = workflow.get_raw
        self.assertEqual(len(raw['x']), len(sx))
        self.assertEqual(len(raw['y']), len(sy))
        self.assertEqual(len(raw['e']), len(se))
        for k in range(len(sx)):
            self.assertEqual(raw['x'][k], sx[k])
            self.assertEqual(raw['y'][k], sy[k])
            self.assertEqual(raw['e'][k], se[k])

    def test_one(self):
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

        # not from Mantid
        self.assertAlmostEqual(results['N1:loglikelihood'][0], -389.94, 2)

        # from Mantid, if not then as a comment
        self.assertAlmostEqual(results['N1:f2.f2.FWHM'][0], 0.055, 3)
        self.assertAlmostEqual(results['N1:f2.f2.beta'][0], 0.794, 3)  # 0.752

        # dont compare amp to Mantid due to different scaling etc.
        self.assertAlmostEqual(results['N1:f2.f2.Amplitude'][0], 0.167, 2)

    def test_two(self):
        """
        Want to check that two calls to the function will append the results
        correctly. So if we use the same input data as above, we expect
        both values to be the same for every item in the dict.
        """
        sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'qse_res.npy'),
                             allow_pickle=True)

        sample = {'x': sx, 'y': sy, 'e': se}
        resolution = {'x': rx, 'y': ry}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, fit_e) = qse_data_main(sample, resolution,
                                             "linear", -0.4, 0.4,
                                             True, results, errors)

        # use the previous results to make it faster
        lbg = LinearBG()
        qse = QSEFunction(lbg, True, rx, ry, -.4, 0.4)
        qse.add_single_SE()
        params = qse.read_from_report(results, 1, 0)

        # call it again
        (results, errors,
         new_x, fit, fit_e) = qse_data_main(sample, resolution,
                                            "linear", -0.4, 0.4,
                                            True, results, errors, params)

        for key in results.keys():
            self.assertEqual(len(results[key]), 2)
            tmp = results[key]
            self.assertAlmostEqual(tmp[0], tmp[1], 3)


if __name__ == '__main__':
    unittest.main()
