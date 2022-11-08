from quasielasticbayes.v2.functions.base import BaseFitFunction
from scipy.optimize import curve_fit
from numpy import ndarray
import numpy as np
from typing import List


def scipy_curve_fit(x_data: ndarray, y_data: ndarray, e_data: ndarray,
                    func: BaseFitFunction,
                    guess: List[float], lower: List[float],
                    upper: List[float]) -> (float, float,
                                            List[float], ndarray):
    """
    A wrapper for the scipy curve fit optimizer.
    :param x_data: x data to fit
    :param y_data: y data to fit
    :param e_data: error data for fit
    :param func: the fitting function
    :param guess: the initial fit values
    :param lower: the lower bounds for the fit params
    :param upper: the upper bounds for the fit params
    :return chi squared, log_10 of the determinant of the Hessian,
    parameters and the y fit values
    """
    params, covar = curve_fit(func, x_data, y_data, guess,
                              sigma=e_data, absolute_sigma=True,
                              maxfev=220000,
                              bounds=(lower, upper))
    fit = func(x_data, *params)

    hessian = np.linalg.inv(covar)
    log_hess_det = np.log10(np.linalg.det(hessian))

    chisq = np.sum(((y_data - fit)/e_data)**2)
    chisq /= (len(x_data) - len(params))

    return chisq, log_hess_det, params, fit
