from quickBayes.fitting.fit_utils import log10_hessian_det
from numpy import ndarray
import numpy as np
from math import exp, log10


def loglikelihood(x_len: ndarray, chi2: float, covar: ndarray,
                  N_peaks: int, beta: float) -> float:
    """
    Calculate the unnormalised logliklihood. log_10 probabilities -> priors

    :param x_len: the length of the x data
    :param chi2: the chi squared for the fit
    :param covar: the covariance matrix
    :param N_peaks: the number of peaks
    :param beta: the scale factor, A_max*(x_max - x_min) eq. 4.17
    :return the loglikelihood
    """

    log_hess_det = log10_hessian_det(covar)
    """
    If the fit is overparameterised the loglikelihood
    will artifically be low. This is because the equation
    has used a Taylor expansion of exp(-0.5*chi^2).
    This requires us to be at a minimum, which is not the
    case with overparameterised data (consider fitting 2 flat
    backgrounds to flat data, there will be a vally of good fits
    that sum to the correct background value).
    We will therefore add a penality to the loglikelihood
    via the hessian.
    """
    if np.max(np.abs(covar)) > 1:
        log_hess_det = 100*np.abs(log_hess_det)

    # want the unscaled chi^2 -> multiple by length of data
    log_chi2 = log10(exp(1.))*chi2*float(x_len)/2.

    log_likelihood = np.sum(
        np.log10(np.asarray([k + 1 for k in range(N_peaks)])))
    log_likelihood += float(N_peaks)*log10(4.*np.pi)
    log_likelihood -= log_chi2
    log_likelihood -= float(N_peaks)*log10(beta)
    log_likelihood -= 0.5*(log_hess_det)

    return log_likelihood
