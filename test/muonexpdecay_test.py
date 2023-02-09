import unittest
from quasielasticbayes.v2.MuonExpDecay import muon_expdecay_main
import numpy as np
import os.path

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
DATA_DIR = os.path.join(DATA_DIR, 'muon')


class MuonExpDecayTest(unittest.TestCase):
    def test_one_decay(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_1.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "linear",
                                                     0.15, 15.0,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -113.05, 2)
        # this is due to a bug ..
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -104.11, 2)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -120.92, 2)
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -130.71, 2)

        self.assertAlmostEqual(results['N1:f2.Amplitude'][0], 0.1, 2)
        self.assertAlmostEqual(results['N1:f2.lambda'][0], 1.04, 2)
        self.assertAlmostEqual(errors['N1:f2.Amplitude'][0], 0.001, 2)
        self.assertAlmostEqual(errors['N1:f2.lambda'][0], 0.032, 2)

    def test_two_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_2.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "linear",
                                                     0.15, 15.0,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -240.91, 2)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -130.99, 2)
        # this is due to a bug ..
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -131.77, 2)
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -129.46, 2)

        self.assertAlmostEqual(results['N2:f2.Amplitude'][0], 0.27, 2)
        self.assertAlmostEqual(results['N2:f2.lambda'][0], 14.58, 2)
        self.assertAlmostEqual(errors['N2:f2.Amplitude'][0], 0.037, 2)
        self.assertAlmostEqual(errors['N2:f2.lambda'][0], 1.25, 2)

        self.assertAlmostEqual(results['N2:f3.Amplitude'][0], 0.1, 2)
        self.assertAlmostEqual(results['N2:f3.lambda'][0], 1.0, 2)
        self.assertAlmostEqual(errors['N2:f3.Amplitude'][0], 0.001, 2)
        self.assertAlmostEqual(errors['N2:f3.lambda'][0], 0.04, 2)

    def test_three_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_3.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "linear",
                                                     0.15, 15.0,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -334.88, 2)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -126.24, 2)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -124.81, 2)
        # this is due to a bug ..
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -123.43, 2)

        self.assertAlmostEqual(results['N3:f2.Amplitude'][0], 0.35, 2)
        self.assertAlmostEqual(results['N3:f2.lambda'][0], 12.23, 2)
        self.assertAlmostEqual(errors['N3:f2.Amplitude'][0], 0.028, 2)
        self.assertAlmostEqual(errors['N3:f2.lambda'][0], 2.06, 2)

        self.assertAlmostEqual(results['N3:f3.Amplitude'][0], 0.09, 2)
        self.assertAlmostEqual(results['N3:f3.lambda'][0], 3.34, 2)
        self.assertAlmostEqual(errors['N3:f3.Amplitude'][0], 0.03, 2)
        self.assertAlmostEqual(errors['N3:f3.lambda'][0], 1.60, 2)

        self.assertAlmostEqual(results['N3:f4.Amplitude'][0], 0.08, 2)
        self.assertAlmostEqual(results['N3:f4.lambda'][0], 0.92, 2)
        self.assertAlmostEqual(errors['N3:f4.Amplitude'][0], 0.02, 2)
        self.assertAlmostEqual(errors['N3:f4.lambda'][0], 0.2, 2)


if __name__ == '__main__':
    unittest.main()
