import unittest
from quasielasticbayes.v2.log_likelihood import loglikelihood


class LoglikelihoodTest(unittest.TestCase):
    def test_1_peak(self):
        x_size = 2231
        log_hess_det = 38.852
        chi2 = 1.314
        N_peaks = 1
        beta = 0.6

        result = loglikelihood(x_size, chi2, log_hess_det, N_peaks, beta)
        self.assertAlmostEqual(result, -654.679, 3)

    def test_2_peaks(self):
        x_size = 2231
        log_hess_det = 36.496
        chi2 = 0.667
        N_peaks = 2
        beta = 0.6

        result = loglikelihood(x_size, chi2, log_hess_det, N_peaks, beta)
        self.assertAlmostEqual(result, -338.437, 3)

    def test_3_peaks(self):
        x_size = 2231
        log_hess_det = 36.253
        chi2 = 0.666
        N_peaks = 3
        beta = 0.6

        result = loglikelihood(x_size, chi2, log_hess_det, N_peaks, beta)
        self.assertAlmostEqual(result, -336.033, 3)


if __name__ == '__main__':
    unittest.main()
