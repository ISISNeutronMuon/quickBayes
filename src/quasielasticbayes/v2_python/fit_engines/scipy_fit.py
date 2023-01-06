from scipy.optimize import curve_fit
from numpy import ndarray
from typing import List, Callable


def scipy_curve_fit(x_data: ndarray, y_data: ndarray, e_data: ndarray,
                    func: Callable,
                    guess: List[float], lower: List[float],
                    upper: List[float]) -> (ndarray, ndarray, ndarray):
    """
    A wrapper for the scipy curve fit optimizer.
    :param x_data: x data to fit
    :param y_data: y data to fit
    :param e_data: error data for fit
    :param func: the fitting function
    :param guess: the initial fit values
    :param lower: the lower bounds for the fit params
    :param upper: the upper bounds for the fit params
    :return parameters, covarience matrix and the y fit values
    """
    max_iterations = 220000
    params, covar = curve_fit(func, x_data, y_data, guess,
                              sigma=e_data, absolute_sigma=True,
                              maxfev=max_iterations,
                              bounds=(lower, upper))

    fit = func(x_data, *params)
    return params, covar, fit
