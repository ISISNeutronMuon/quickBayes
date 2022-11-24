from numpy import ndarray
import numpy as np
from typing import Callable
from scipy.stats import t


def log10_hessian_det(covar: ndarray) -> float:
    """

    """
    hessian = np.linalg.inv(covar)
    return np.log10(np.linalg.det(hessian))


def chi_squared(x_data: ndarray, y_data: ndarray, e_data: ndarray,
                fit: ndarray, params: ndarray) -> float:
    """
    Method for calculating the reduced chi squared value
    :param x_data: x data to fit
    :param y_data: y data to fit
    :param e_data: error data for fit
    :param fit: the fit
    :param params: the fit parameters
    :return the chi squared value
    """
    chisq = np.sum(((y_data - fit)/e_data)**2)
    chisq /= (len(x_data) - len(params))
    return chisq


def param_errors(covar: ndarray) -> ndarray:
    """
    Get the errors for the parameters
    :param covar: the covarience matrix
    :return the errors for the parameters
    """
    return np.sqrt(np.diag(covar))


def derivative(x_data: ndarray, params: ndarray, func: Callable) -> ndarray:
    """
    Get numerical derivative for a function
    :param x_data: the x data
    :param params: the paramaters
    :param func: the function
    :return numerical derivatives (with respect to fitting parameter)
    """
    df_by_dp = []
    N = len(params)
    for j in range(N):
        # only want to change one parameter at a time
        dparams = np.zeors(N)
        # small (0.1%) change in parameter value
        dparams[j] = params[j]*0.001
        # forward difference
        df_by_dp.append((func(x, *(params + dparams)) -
                         func(x, *(params)))/np.sum(dparams))
    return df_by_dp


def fit_errors(x_data: ndarray, params: ndarray, fit: ndarray,
               covar: ndarray, df_by_dp: ndarray) -> ndarray:
    """

    """
    confidence = 0.6826  # 2 sigma

    prob = 0.5 + confidence/2.  # even distribution above and below data point
    N = len(params)
    M = len(x_data)
    dof = M - N
    tval = t.ppf(prob, dof)

    df_sq = np.zeros(M)
    for j in range(N):
        for k in range(N):
            df_sq += df_by_dp[j]*df_by_dp[k]*covar[j, k]
    df = np.sqrt(df_sq)

    return tval*df
