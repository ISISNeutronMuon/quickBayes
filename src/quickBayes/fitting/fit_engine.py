from numpy import ndarray
from typing import Callable
from abc import abstractmethod
import numpy as np
from quickBayes.fitting.fit_utils import (chi_squared,
                                          param_errors,
                                          derivative,
                                          fit_errors,
                                          var, res)
from quickBayes.utils.spline import spline


class FitEngine(object):
    """
    A basic class for the fit engine, includes a history
    """
    def __init__(self, name: str, x_data: ndarray, y_data: ndarray,
                 e_data: ndarray):
        """
        Creates the basic fit engine class
        Stores useful information about each fit
        :param name: name of the fit engine
        :param x_data: original x data (can fit to an interpolation)
        :param y_data: original y data (can fit to an interpolation)
        :param e_data: original e data (can fit to an interpolation)
        """
        self._name = name
        self._fits = []
        self._fit_errors = []
        self._params = []
        self._param_errors = []
        self._msgs = []
        self._chi2 = []
        self._covars = []
        self._diffs = []
        self._fit = None

        self._x_data = x_data
        self._y_data = y_data
        self._e_data = e_data

    def get_chi_squared(self, index: int = -1) -> float:
        """
        Get the chi squared value
        :param index: the index (number) of chi squared that you want,
        count from 0
        :return chi squared value
        """
        return self._chi2[index]

    def get_covariance_matrix(self, index: int = -1) -> ndarray:
        """
        :param index: the index (number) of covariance matrix that you want,
        count from 0
        :return covariance matrix
        """
        return self._covars[index]

    def get_fit_values(self, index: int = -1) -> (ndarray, ndarray,
                                                  ndarray, ndarray, ndarray):
        """
        Gets the fit values from the history.
        These are interpolated onto the same x values as the
        fit history was created with.
        The fit errors are calculated by quadrature sqrt(sigma_f^2 + sigma_y^2)
        :param index: the index (number) of fit that you want,
        counts from 0
        :return fit values (x data, y values, y errors, diffs, diff errors)
        """
        return (self._x_data, self._fits[index], self._fit_errors[index],
                self._diffs[index], np.sqrt(self._fit_errors[index]**2 +
                                            self._e_data**2))

    def get_fit_parameters(self, index: int = -1) -> (ndarray, ndarray):
        """
        Get the fit parameters and errors, that you want
        :param index: the index (number) of fit parameters (and errors) that
        you want, counts from 0
        :return list of fit parameters and their errors
        """
        return self._params[index], self._param_errors[index]

    @property
    def name(self) -> str:
        """
        :return name of the fit engine
        """
        return self._name

    @abstractmethod
    def _do_fit(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
                func: Callable) -> ndarray:
        """
        The main function for doing the fitting_algorithm
        :param x_data: the x data to fit
        :param y_data: the y data to fit
        :param e_data: the error data to fit
        :return the fit parameters
        """
        raise NotImplementedError()

    def add_fit(self, x_data: ndarray, func: Callable, df_by_dp: ndarray,
                params: ndarray) -> None:
        """
        Adds the fit result to the fit history
        :param x_data: the x data to fit
        :param func: the fitting function
        :param df_by_dp: the derivatives wrt the parameters
        :param params: the parameters from the fit
        """
        fit_y = func(x_data, *params)
        errors = fit_errors(x_data, params, fit_y, self._covars[-1], df_by_dp)
        self._fit = fit_y
        # record fit on same x axis as the FitHistory was create with
        if not np.array_equal(x_data, self._x_data):
            fit_y = spline(x_data, fit_y, self._x_data)
            errors = spline(x_data, errors, self._x_data)
        self._fits.append(fit_y)
        self._fit_errors.append(errors)
        self._diffs.append(fit_y - self._y_data)

    def add_params(self, params) -> None:
        """
        Add the parameters and errors to the history
        :param params: the fit parameters
        """
        self._params.append(params)
        self._param_errors.append(param_errors(self._covars[-1]))

    def do_fit(self, x_data: ndarray, y_data: ndarray, e_data: ndarray,
               func: Callable) -> None:
        """
        Call for doing a fit and updating the history
        :param x_data: the x data to fit against
        :param y_data: the y data to fit against
        :param e_data: the error data to fit against
        :param func: the fitting function
        """
        params = self._do_fit(x_data, y_data, e_data, func)

        df_by_dp = derivative(x_data, params, func)
        self.calculate_covar(x_data, y_data, e_data, func, df_by_dp, params)
        self.add_params(params)
        self.add_fit(x_data, func, df_by_dp, params)
        self._chi2.append(chi_squared(x_data, y_data, e_data,
                                      self._fit, params))
        self._fit = None

    def calculate_covar(self, x_data: ndarray, y_data: ndarray,
                        e_data: ndarray,
                        func: Callable, df_by_dp: ndarray,
                        params: ndarray) -> None:
        """
        Calculate the covariance matrix and add it to the history
        :param x_data: the x data to fitted against
        :param y_data: the y data to fitted against
        :param e_data: the error data to fitted against
        :param func: the fitting function
        :param df_by_dp: the derivatives wrt the parameters
        :param params: the fit parameters
        """
        # make Jacobian matrix (derivatives)
        jac = np.array(df_by_dp).T
        # factorize the matrix
        _, upper_triangle = np.linalg.qr(jac)
        # Calculate the inverse value of upper triangle
        inverse = np.linalg.solve(upper_triangle,
                                  np.identity(upper_triangle.shape[0]))
        # Matrix multiplication: (J^T J)^{-1}
        JTJ_inv = np.matmul(inverse, inverse.transpose())

        # weight: sum( y - f)^2/sum( (y-f)^2/e^2) -> cannot cancel due to sum
        weight = var(func, x_data, y_data, params)/res(func, x_data, y_data,
                                                       e_data, params)
        CovMatrix = JTJ_inv * weight
        self._covars.append(CovMatrix)
