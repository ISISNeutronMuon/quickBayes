from gofit import multistart
from numpy import ndarray
from typing import Callable
from quickBayes.fitting.fit_engine import FitEngine


"""
This file contains all of the code needed for a gofit
engine, because of the way gofit works we need to provide
a cost function. However, when the cost function is called
by gofit it is not provided with the original data. Hence,
we need to construct a cost function class first that has
the data as a member variable.
"""


class ChiSquared(object):
    def __init__(self, x_data: ndarray, y_data: ndarray,
                 e_data: ndarray, func: Callable):
        """
        A chi^2 cost function class for use with gofit
        :param x_data: x data that fitted against
        :param y_data: y data that fitted against
        :param e_data: e data that fitted against
        :param func: the fitting function used
        """
        self._x_data = x_data
        self._y_data = y_data
        self._e_data = e_data
        self._func = func

    def __call__(self, params: ndarray) -> float:
        """
        Calls the evaluation of the cost function
        :param params: the fit parameters
        :return the cost function evaluation
        """
        fit = self._func(self._x_data, *params)
        return (fit - self._y_data)**2 / self._e_data**2


class GoFitEngine(FitEngine):
    """
    A gofit multistart fit engine.
    This will use gofit's multistart to
    fit data.
    """

    def __init__(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                 lower: ndarray, upper: ndarray, samples: int = 10,
                 max_iterations: int = 220000):
        """
        Creates the scipy curve fit engine class
        Stores useful information about each fit
        :param name: name of the fit engine
        :param x_data: original x data (can fit to an interpolation)
        :param y_data: original y data (can fit to an interpolation)
        :param e_data: original e data (can fit to an interpolation)
        :param lower: the lower bounds for the fit parameters
        :param upper: the upper bounds for the fit parameters
        :param samples: the number of samples to use in multistart
        :param max_iterations: the maximum number of iterations for the fit
        """
        super().__init__("gofit", x_data, y_data, e_data)
        # extra parameters
        self.set_bounds_and_N_params(lower, upper)
        self._max_iterations = max_iterations
        self._samples = samples

    def set_bounds_and_N_params(self, lower: ndarray, upper: ndarray) -> None:
        """
        Sets the current bounds and number of parameters for the fit function.
        If the functional form changes this method will need to be called
        with updated values.
        :param lower: the lower bound for the function parameters
        :param upper: the upper bound for the function parameters
        """
        # validate
        if len(upper) != len(lower):
            raise ValueError(f"The lower {lower} and "
                             f"upper {upper} bounds must "
                             "be the same length")
        self._lower = lower
        self._upper = upper
        self._N_params = len(upper)

    def _do_fit(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                func: Callable) -> ndarray:
        """
        Calls gofit multistart
        :param x_data: the x data to fit
        :param y_data: the y data to fit
        :param e_data: the error data to fit
        :param func: the fitting function
        :return the fit parameters
        """
        cost_function = ChiSquared(x_data, y_data, e_data, func)

        data_length = len(x_data)

        params, _ = multistart(data_length, self._N_params,
                               self._lower, self._upper,
                               cost_function, samples=self._samples,
                               maxit=self._max_iterations)

        return params
