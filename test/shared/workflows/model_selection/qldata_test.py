import unittest
from quickBayes.workflow.model_selection.QlData import ql_data_main, QLData
from quickBayes.functions.qldata_function import QlDataFunction
from quickBayes.functions.BG import LinearBG
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')


"""
Paramater result are from Mantid v6.5 on Windows
"""


class QlDataTest(unittest.TestCase):

    def test_raw(self):
        sx = [0, 1, 3, 5, 6.9]
        sy = [1, 2, .2, .4, .1]
        se = [.1, .2, .3, .4, .5]

        sample = {'x': sx, 'y': sy, 'e': se}
        res = {'x': sx, 'y': sy}

        results = {}
        errors = {}
        workflow = QLData(results, errors)
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
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'resolution_data_red.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        resolution = {'x': rx, 'y': ry}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = ql_data_main(sample, resolution,
                                               "linear", -0.4, 0.4,
                                               True, results, errors)

        # not from Mantid
        self.assertAlmostEqual(results['N1:loglikelihood'][0], -659.4, 2)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -349.17, 2)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -350.87, 2)

        # from Mantid, if not then as a comment
        self.assertAlmostEqual(results['N1:f2.f2.Gamma'][0], 0.0566, 2)

        self.assertAlmostEqual(results['N2:f2.f2.Gamma'][0], 0.201, 1)  # 0.326
        self.assertAlmostEqual(results['N2:f2.f3.Gamma'][0], 0.04, 2)   # 0.046

        self.assertAlmostEqual(results['N3:f2.f2.Gamma'][0], 0.216, 2)  # 0.347
        self.assertAlmostEqual(results['N3:f2.f3.Gamma'][0], 0.050, 1)
        self.assertAlmostEqual(results['N3:f2.f4.Gamma'][0], 0.016, 2)  # 0.012

        self.assertAlmostEqual(errors['N1:f2.f2.Gamma'][0], 0.00018, 5)

        self.assertAlmostEqual(errors['N2:f2.f2.Gamma'][0], 0.010, 2)
        self.assertAlmostEqual(errors['N2:f2.f3.Gamma'][0], 0.0006, 4)

        self.assertAlmostEqual(errors['N3:f2.f2.Gamma'][0], 0.02, 2)
        self.assertAlmostEqual(errors['N3:f2.f3.Gamma'][0], 0.006, 3)
        self.assertAlmostEqual(errors['N3:f2.f4.Gamma'][0], 0.03, 2)

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

        self.assertAlmostEqual(errors['N1:f2.f2.EISF'][0], 0.006, 3)

        self.assertAlmostEqual(errors['N2:f2.f2.EISF'][0], 0.02, 2)
        self.assertAlmostEqual(errors['N2:f2.f3.EISF'][0], 0.011, 3)

        self.assertAlmostEqual(errors['N3:f2.f2.EISF'][0], 0.4, 1)
        self.assertAlmostEqual(errors['N3:f2.f3.EISF'][0], 0.116, 2)
        self.assertAlmostEqual(errors['N3:f2.f4.EISF'][0], 0.03, 2)

    def test_two(self):
        """
        Want to check that two calls to the function will append the results
        correctly. So if we use the same input data as above, we expect
        both values to be the same for every item in the dict.
        """
        sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
        rx, ry, re = np.load(os.path.join(DATA_DIR, 'resolution_data_red.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        resolution = {'x': rx, 'y': ry}
        results = {}
        errors = {}

        (results, errors,
         new_x, fit, fit_errors) = ql_data_main(sample, resolution,
                                                "linear", -0.4, 0.4,
                                                True, results, errors)

        # call it again
        ql = QlDataFunction(LinearBG(), True, rx, ry, -0.4, 0.4)
        ql.add_single_lorentzian()
        params = ql.read_from_report(results, 1, -1)

        (results, errors,
         new_x, fit2, fit_errors2) = ql_data_main(sample, resolution,
                                                  "linear", -0.4, 0.4,
                                                  True, results,
                                                  errors, params)

        params = ql.read_from_report(results, 1, -1)
        for key in results.keys():
            self.assertEqual(len(results[key]), 2)
            tmp = results[key]

            percentage_change = 100.*np.abs((tmp[0] - tmp[1])/tmp[0])
            self.assertLessEqual(percentage_change, 20.)

        for key in errors.keys():
            self.assertEqual(len(results[key]), 2)
            tmp = errors[key]
            self.assertAlmostEqual(tmp[0], tmp[1], 3)

            percentage_change = 100.*np.abs((tmp[0] - tmp[1])/tmp[0])
            self.assertLessEqual(percentage_change, 20.)

        for j in range(3):
            for k in range(len(fit[j])):
                self.assertAlmostEqual(fit[j][k], fit2[j][k], 3)
                self.assertAlmostEqual(fit_errors[j][k], fit_errors2[j][k], 3)


if __name__ == '__main__':
    unittest.main()
