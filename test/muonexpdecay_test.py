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
                                                     0.2, 15.,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -104., 0)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -463., 0)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -1857., 0)
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -547., 0)

        self.assertAlmostEqual(results['N1:f2.Amplitude'][0], 0.1, 2)
        self.assertAlmostEqual(results['N1:f2.lambda'][0], 1.03, 2)
        self.assertAlmostEqual(errors['N1:f2.Amplitude'][0], 0.001, 2)
        self.assertAlmostEqual(errors['N1:f2.lambda'][0], 0.032, 2)

    def test_two_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_2.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "linear",
                                                     0.16, 14.5,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -142, 0)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -114., 0)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -712, 0)
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -298, 0)

        self.assertAlmostEqual(results['N2:f2.Amplitude'][0], 0.17, 2)
        self.assertAlmostEqual(results['N2:f2.lambda'][0], 3.87, 2)
        self.assertAlmostEqual(errors['N2:f2.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f2.lambda'][0], 0.39, 2)

        self.assertAlmostEqual(results['N2:f3.Amplitude'][0], 0.1, 2)
        self.assertAlmostEqual(results['N2:f3.lambda'][0], 1.02, 2)
        self.assertAlmostEqual(errors['N2:f3.Amplitude'][0], 0.01, 2)
        self.assertAlmostEqual(errors['N2:f3.lambda'][0], 0.12, 2)

    def test_three_decays(self):
        sx, sy, se = np.loadtxt(os.path.join(DATA_DIR, 'muon_expdecay_3.npy'))

        sample = {'x': sx, 'y': sy, 'e': se}
        results = {}
        errors = {}

        (results, errors,
         new_x, fits, f_errors) = muon_expdecay_main(sample, "linear",
                                                     0.15, 15.0,
                                                     results, errors)

        self.assertAlmostEqual(results['N1:loglikelihood'][0], -201, 0)
        self.assertAlmostEqual(results['N2:loglikelihood'][0], -118, 0)
        self.assertAlmostEqual(results['N3:loglikelihood'][0], -1569, 0)
        self.assertAlmostEqual(results['N4:loglikelihood'][0], -1469, 0)

        self.assertAlmostEqual(results['N3:f2.Amplitude'][0], 0.45, 2)
        self.assertAlmostEqual(results['N3:f2.lambda'][0], 13.78, 2)
        self.assertAlmostEqual(errors['N3:f2.Amplitude'][0], 0.20, 2)
        self.assertAlmostEqual(errors['N3:f2.lambda'][0], 4.78, 1)

        self.assertAlmostEqual(results['N3:f3.Amplitude'][0], 0.09, 2)
        self.assertAlmostEqual(results['N3:f3.lambda'][0], 3.51, 2)
        self.assertAlmostEqual(errors['N3:f3.Amplitude'][0], 0.04, 2)
        self.assertAlmostEqual(errors['N3:f3.lambda'][0], 1.90, 2)

        self.assertAlmostEqual(results['N3:f4.Amplitude'][0], 0.08, 2)
        self.assertAlmostEqual(results['N3:f4.lambda'][0], 0.95, 2)
        self.assertAlmostEqual(errors['N3:f4.Amplitude'][0], 0.02, 2)
        self.assertAlmostEqual(errors['N3:f4.lambda'][0], 0.2, 2)


if __name__ == '__main__':
    unittest.main()
