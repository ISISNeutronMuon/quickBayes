from quickBayes.fitting.fit_utils import log10_hessian_det
from numpy import ndarray
import numpy as np
from math import exp, log10


def loglikelihood(x_len: ndarray, chi2: float, covar: ndarray,
                  N_peaks: int, beta: float) -> float:
    """
    Calculate the unnormalised logliklihood.
    log_10 probabilities -> priors
    The equation for the probability is taken from
    "Data Analysis A bayesian tutorial second edition",
    by D. S. Sivia
    equation 4.20 page 88:
    P(M|{D_k},I) prop \\frac{N! (4\\pi)^N}{\\beta^N\\sqrt(\\det(H))}\\exp{(-\\frac{\\chi^2}{2})}  # noqa E501
    where \\beta = (x_{max}-x_{min})A_{max}, H is the Hessian matrix
    We want this as a logliklihood, so take \\log_{10} and use:
    - \\log_{10}(\\exp{\\alpha}) = \\ln(\\exp{\\alpha})\\log_{10}(\\exp{1})
    \\log_{10}(\\exp({\\alpha}) = \\alpha\\log_{10}(\\exp{1})
    - N! = \\Pi_{i=0}{N} i => \\log{N!} = \\sum_{i=0}^N \\log(i)
    - \\log(ab) = \\log(a) + \\log(b)
    - \\log(\\frac{a}{b}) = \\log(a) - \\log(b)
    - \\log(a^N) = N\\log(a)
    to get:
    \\sum_{j=1}{N}\\log(j) + N\\log(4\\pi) - 0.5\\chi^2\\log(exp(1))
    - N\\log(\\beta) - 0.5\\log(\\det(H))
    :param x_len: the length of the x data
    :param chi2: the chi squared for the fit
    :param covar: the covariance matrix
    :param N_peaks: the number of peaks
    :param beta: the scale factor, A_max*(x_max - x_min) eq. 4.17
    :return the loglikelihood
    """

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
    log_hess_det = log10_hessian_det(covar)
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
