from numpy import ndarray
import numpy as np
from typing import Callable
from scipy.stats import t as student_t_dist


TWO_SIGMA = 0.6826


def var(func: Callable, x_data: ndarray, y_data: ndarray,
        params: ndarray) -> float:
    """
    Calculate the variance of the data: sum_i (f(x_i) - y_i)^2
    :param func: the fitting function
    :param x_data: the x data
    :param y_data: the y data
    :param params: the fitting parameters
    :return the variance
    """
    return np.sum(pow(func(x_data, *params) - y_data, 2))


def res(func, x_data, y_data, e_data, param):
    """
    Calculate the residuals of the data: sum_i ((f(x_i) - y_i)/e_i)^2
    :param func: the fitting function
    :param x_data: the x data
    :param y_data: the y data
    :param e_data: the e data
    :param params: the fitting parameters
    :return the residuals
    """
    tmp = (func(x_data, *param) - y_data)/e_data
    return np.sum(tmp*tmp)


def log10_hessian_det(covar: ndarray) -> float:
    """
    Calculate the log base 10 of the determinant
    of the Hessian matrix
    :param covar: the covarience matrix
    :return the log of the determinant of the
    Hessian matrix
    """
    hessian = np.linalg.inv(covar)
    det = np.linalg.det(hessian)
    # cannot have a negative value in log
    if det <= 0:
        det = 1.e-9
    return np.log10(det)


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
        dparams = np.zeros(N)
        # small (0.1%) change in parameter value
        dparams[j] = params[j]*0.001
        # forward difference
        df_by_dp.append((func(x_data, *(params + dparams)) -
                         func(x_data, *(params)))/np.sum(dparams))
    return df_by_dp


def fit_errors(x_data: ndarray, params: ndarray, fit: ndarray,
               covar: ndarray, df_by_dp: ndarray) -> ndarray:
    """
    Generate the errors for the fit line
    :param x_data: the x data
    :param params: the parameters for the function
    :param fit: the y data values from the fit
    :param covar: the covarience matrix
    :param df_by_dp: the derivatives
    :return the error values
    """
    confidence = TWO_SIGMA

    prob = 0.5 + confidence/2.  # even distribution above and below data point
    N = len(params)
    M = len(x_data)
    dof = M - N
    tval = student_t_dist.ppf(prob, dof)

    df_sq = np.zeros(M)
    for j in range(N):
        for k in range(N):
            df_sq += df_by_dp[j]*df_by_dp[k]*covar[j, k]
    df = np.sqrt(df_sq)

    return tval*df
