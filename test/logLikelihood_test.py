import unittest
from unittest import mock
from quickBayes.log_likelihood import loglikelihood


class LoglikelihoodTest(unittest.TestCase):
    @mock.patch("quickBayes.log_likelihood.log10_hessian_det")
    def test_1_peak(self, mock_log_hess):
        x_size = 2231
        mock_log_hess.return_value = 38.852
        chi2 = 1.314
        N_peaks = 1
        beta = 0.6
        covar = [[.1, .2], [.1, -.4]]
        result = loglikelihood(x_size, chi2, covar, N_peaks, beta)
        self.assertAlmostEqual(result, -654.679, 3)

    @mock.patch("quickBayes.log_likelihood.log10_hessian_det")
    def test_2_peaks(self, mock_log_hess):
        x_size = 2231
        mock_log_hess.return_value = 36.496
        chi2 = 0.667
        N_peaks = 2
        beta = 0.6
        covar = [[.1, .2], [.1, -.4]]

        result = loglikelihood(x_size, chi2, covar, N_peaks, beta)
        self.assertAlmostEqual(result, -338.437, 3)

    @mock.patch("quickBayes.log_likelihood.log10_hessian_det")
    def test_3_peaks(self, mock_log_hess):
        x_size = 2231
        mock_log_hess.return_value = 36.253
        chi2 = 0.666
        N_peaks = 3
        beta = 0.6
        covar = [[.1, .2], [.1, -.4]]

        result = loglikelihood(x_size, chi2, covar, N_peaks, beta)
        self.assertAlmostEqual(result, -336.033, 3)

    @mock.patch("quickBayes.log_likelihood.log10_hessian_det")
    def test_over_optimized(self, mock_log_hess):
        x_size = 2231
        mock_log_hess.return_value = 36.253
        chi2 = 0.666
        N_peaks = 3
        beta = 0.6
        covar = [[.1, 1.1], [1.1, 0.4]]

        result = loglikelihood(x_size, chi2, covar, N_peaks, beta)
        self.assertAlmostEqual(result, -2131, 0)

    @mock.patch("quickBayes.log_likelihood.log10_hessian_det")
    def test_over_optimized_negative(self, mock_log_hess):
        x_size = 2231
        mock_log_hess.return_value = 36.253
        chi2 = 0.666
        N_peaks = 3
        beta = 0.6
        covar = [[.1, 0.1], [0.1, -1.4]]

        result = loglikelihood(x_size, chi2, covar, N_peaks, beta)
        self.assertAlmostEqual(result, -2131, 0)


if __name__ == '__main__':
    unittest.main()
