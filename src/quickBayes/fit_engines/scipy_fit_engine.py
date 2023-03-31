from scipy.optimize import curve_fit
from numpy import ndarray
from typing import Callable
from quickBayes.fitting.fit_engine import FitEngine


class ScipyFitEngine(FitEngine):
    """
    A scipy curve fit fit engine.
    This will use scipy's curve fit to
    fit data.
    """

    def __init__(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                 lower: ndarray, upper: ndarray, guess: ndarray,
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
        :param guess: the initial guess for the fit parameters
        :param max_iterations: the maximum number of iterations for the fit
        """
        super().__init__("scipy", x_data, y_data, e_data)
        # extra parameters
        self.set_guess_and_bounds(guess, lower, upper)
        self._max_iterations = max_iterations

    def set_guess_and_bounds(self, guess: ndarray,
                             lower: ndarray, upper: ndarray) -> None:
        """
        Sets the current guess and bounds for the fit function.
        If the functional form changes this method will need to be called
        with updated values.
        :param guess: the initial guess for the fit function parameters
        :param lower: the lower bound for the function parameters
        :param upper: the upper bound for the function parameters
        """
        # validate
        if len(guess) != len(upper) or len(upper) != len(lower):
            raise ValueError(f"The guess {guess}, lower {lower} and "
                             f"upper {upper} bounds must "
                             "be the same length")
        self._guess = guess
        self._lower = lower
        self._upper = upper

    def _do_fit(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                func: Callable) -> ndarray:
        """
        Calls scipy curve fit
        :param x_data: the x data to fit
        :param y_data: the y data to fit
        :param e_data: the error data to fit
        :return the fit parameters
        """
        params, covar = curve_fit(func, x_data, y_data, self._guess,
                                  sigma=e_data, absolute_sigma=True,
                                  maxfev=self._max_iterations,
                                  bounds=(self._lower, self._upper))
        self._covars.append(covar)
        return params

    def calculate_covar(self, x_data: ndarray, y_data: ndarray,
                        e_data: ndarray,
                        func: Callable, df_by_dp: ndarray,
                        params: ndarray) -> None:
        """
        Calculated the covariance matrix during the fit
        So do nothing here as we already have the value
        :param x_data: the x data to fitted against
        :param y_data: the y data to fitted against
        :param e_data: the error data to fitted against
        :param func: the fitting function
        :param df_by_dp: the derivatives wrt the parameters
        :param params: the fit parameters
        """
        return
