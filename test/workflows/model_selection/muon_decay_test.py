import unittest
from quickBayes.workflow.model_selection.muon_decay import muon_expdecay_main
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data')
DATA_DIR = os.path.join(DATA_DIR, 'muon')


class MuonExpDecayTest(unittest.TestCase):
    def test_one_decay(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_1.npy'))
        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "flat",
                                                     0.12, 15.,
                                                     results, errors)
        # just need to check that the correct answer is the most likely
        # so all other options should be more negative
        expected = results['N1:loglikelihood'][0]
        self.assertAlmostEqual(expected, -100., 0)
        self.assertLess(results['N2:loglikelihood'][0], expected)
        self.assertLess(results['N3:loglikelihood'][0], expected)
        self.assertLess(results['N4:loglikelihood'][0], expected)

        self.assertAlmostEqual(results['N1:f2.Amplitude'][0], 0.1, 2)
        self.assertAlmostEqual(results['N1:f2.lambda'][0], 1.03, 2)
        self.assertAlmostEqual(errors['N1:f2.Amplitude'][0], 0.001, 3)
        self.assertAlmostEqual(errors['N1:f2.lambda'][0], 0.02, 2)

    def test_two_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_2.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "flat",
                                                     0.16, 15,
                                                     results, errors)

        expected = results['N2:loglikelihood'][0]
        self.assertLess(results['N1:loglikelihood'][0], expected)
        self.assertAlmostEqual(expected, -114., 0)
        self.assertLess(results['N3:loglikelihood'][0], expected)
        self.assertLess(results['N4:loglikelihood'][0], expected)

        self.assertAlmostEqual(results['N2:f2.Amplitude'][0], 0.17, 2)
        self.assertAlmostEqual(results['N2:f2.lambda'][0], 3.86, 2)
        self.assertAlmostEqual(errors['N2:f2.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f2.lambda'][0], 0.34, 2)

        self.assertAlmostEqual(results['N2:f3.Amplitude'][0], 0.10, 2)
        self.assertAlmostEqual(results['N2:f3.lambda'][0], 1.01, 2)
        self.assertAlmostEqual(errors['N2:f3.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f3.lambda'][0], 0.08, 2)

    def test_three_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_3.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        """
        This data was generated with 3 decays. However,
        since this method gives the loglikelihood given
        the data, it is possible for less decays to be more likely.
        This is due to noise making it difficult to identify 3
        destinct decays.
        This used 20 MEvents
        """

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "flat",
                                                     0.2, 14.0,
                                                     results, errors)

        expected = results['N2:loglikelihood'][0]
        self.assertLess(results['N1:loglikelihood'][0], expected)
        self.assertAlmostEqual(expected, -106, 0)
        self.assertLess(results['N3:loglikelihood'][0], expected)
        self.assertLess(results['N4:loglikelihood'][0], expected)

        self.assertAlmostEqual(results['N2:f2.Amplitude'][0], 0.32, 2)
        self.assertAlmostEqual(results['N2:f2.lambda'][0], 5.05, 2)
        self.assertAlmostEqual(errors['N2:f2.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f2.lambda'][0], 0.30, 1)

        self.assertAlmostEqual(results['N2:f3.Amplitude'][0], 0.12, 2)
        self.assertAlmostEqual(results['N2:f3.lambda'][0], 1.13, 2)
        self.assertAlmostEqual(errors['N2:f3.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f3.lambda'][0], 0.07, 2)

    def test_three_decays_big(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR,
                                             'muon_expdecay_3_big.npy'))
        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        """
        This data was generated with 3 decays. However,
        since this method gives the loglikelihood given
        the data, it is possible for less decays to be more likely.
        This is due to noise making it difficult to identify 3
        destinct decays.
        This used 80 MEvents to get good stats, to pick out the
        3 decays
        """

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "flat",
                                                     0.11, 15.0,
                                                     results, errors)
        expected = results['N3:loglikelihood'][0]
        self.assertLess(results['N1:loglikelihood'][0], expected)
        self.assertLess(results['N2:loglikelihood'][0], expected)
        self.assertAlmostEqual(expected, -101, 0)
        self.assertLess(results['N4:loglikelihood'][0], expected)

        self.assertAlmostEqual(results['N3:f2.Amplitude'][0], 0.40, 2)
        self.assertAlmostEqual(results['N3:f2.lambda'][0], 6.17, 2)
        self.assertAlmostEqual(errors['N3:f2.Amplitude'][0], 0.08, 2)
        self.assertAlmostEqual(errors['N3:f2.lambda'][0], 0.62, 1)

        self.assertAlmostEqual(results['N3:f3.Amplitude'][0], 0.24, 2)
        self.assertAlmostEqual(results['N3:f3.lambda'][0], 2.68, 2)
        self.assertAlmostEqual(errors['N3:f3.Amplitude'][0], 0.06, 2)
        self.assertAlmostEqual(errors['N3:f3.lambda'][0], 0.66, 2)

        self.assertAlmostEqual(results['N3:f4.Amplitude'][0], 0.14, 2)
        self.assertAlmostEqual(results['N3:f4.lambda'][0], 0.98, 2)
        self.assertAlmostEqual(errors['N3:f4.Amplitude'][0], 0.03, 2)
        self.assertAlmostEqual(errors['N3:f4.lambda'][0], 0.09, 2)


if __name__ == '__main__':
    unittest.main()
