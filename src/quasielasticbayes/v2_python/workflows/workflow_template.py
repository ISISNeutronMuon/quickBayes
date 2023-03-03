from quasielasticbayes.v2.fitting.scipy_engine import ScipyFitEngine
from quasielasticbayes.v2.functions.base import BaseFitFunction

from quasielasticbayes.v2.utils.general import update_guess
from quasielasticbayes.v2.log_likelihood import loglikelihood

from numpy import ndarray
import numpy as np
from typing import Dict
from abc import abstractmethod


class Workflow(object):
    """
    This is a class for the quick bayes workflow.
    Each method can be overwritten to provide
    unique functionality for the specific
    use case.

    The inherited class must include:
    - _update_function method to increment the fit function

    The properties are:
    - fit_engine
    - get_parameters_and_errors

    To add a fit engine:
    - set_scipy_engine (scipy curve fit)

    Other methods:
    - preprocess_data
    - update_fit_engine
    - update_function (call this one not the overwritten one)
    - report
    - execute
    - update_scipy_fit_engine
    """
    def __init__(self, results: Dict[str, ndarray],
                 results_errors: Dict[str, ndarray]):
        """
        Set the results and error dicts for reporting
        :param results: dict of parameter values
        :param results_errors: dict of parameter errors
        """
        self._engine = None
        self._results_dict = results
        self._errors_dict = results_errors
        self._data = None

    @property
    def fit_engine(self):
        """
        Simple method for getting the fit engine
        :reeturn the fit engine used
        """
        return self._engine

    @property
    def get_parameters_and_errors(self) -> (Dict[str, float],
                                            Dict[str, float]):
        """
        Method to get the dict's of the fit parameters and errors.
        :return dict of fit parameters, dict of fit parameter errors
        """
        return self._results_dict, self._errors_dict

    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        *args: float) -> None:
        """
        The preprocessing needed for the data.
        This simple case just assigns the data values.
        This can be overwritten in a derived class to
        include extra processing (e.g. cropping or
        interpolating the data).
        :param x_data: the x data to fit to
        :param y_data: the y data to fit to
        :param e_data: the errors for the y data
        :param *args: additional arguments that might be needed
        """
        self._data = {'x': x_data, 'y': y_data, 'e': e_data}

    def update_fit_engine(self, func: BaseFitFunction, params: ndarray,
                          *args: float) -> None:
        """
        This updates the fit engine specific properties.
        e.g. the bounds and guess for scipy
        :param func: the fitting function
        :param params: the fitting parameters
        :param *args: additional arguments
        """
        if self._engine.name == 'scipy':
            self.update_scipy_fit_engine(func, params)
        return

    @abstractmethod
    def _update_function(self, func: BaseFitFunction) -> BaseFitFunction:
        """
        This method updates the fitting function when the number of
        features has been incremented.
        It will be unique to the workflow.
        :param func: the fitting function that needs modifing
        :return the modified fitting function
        """
        raise NotImplementedError()

    def update_function(self, func: BaseFitFunction,
                        N: int) -> BaseFitFunction:
        """
        This method updates the fitting function when the number of
        features has been incremented.
        It will be unique to the workflow.
        :param func: the fitting function that needs modifing
        :return the modified fitting function
        """
        function = self._update_function(func)
        function.update_prefix(f'N{N}:')
        return function

    def report(self, func: BaseFitFunction, N: int, beta: float) -> ndarray:
        """
        Reports the latest fit parameters and records the fit
        parameters and their errors into dicts.
        :param func: the fitting function used
        :param N: the number of features used
        :param beta: the beta scaling factor
        :return the fit parameters
        """
        params, errors = self._engine.get_fit_parameters()

        self._results_dict = func.report(self._results_dict, *params)
        self._errors_dict = func.report_errors(self._errors_dict,
                                               errors, params)

        prob_name = f'N{N}:loglikelihood'
        n_data = len(self._data['y'])
        chi2 = self._engine.get_chi_squared()
        covar = self._engine.get_covariance_matrix()

        if prob_name in self._results_dict.keys():
            self._results_dict[prob_name].append(loglikelihood(n_data,
                                                               chi2,
                                                               covar,
                                                               N, beta))
        else:
            self._results_dict[prob_name] = [loglikelihood(n_data,
                                                           chi2,
                                                           covar,
                                                           N, beta)]

        return params

    def execute(self, max_num_features: int,
                func: BaseFitFunction,
                params: ndarray = []) -> BaseFitFunction:
        """
        The main part of the analysis.
        It increments the number of features in the fitting function,
        does a fit and then records the results.
        :param max_num_features: the maximum number of features
        :param func: the fitting function
        :param params: the (optional) initial fit parameters
        :return the update fitting function
        """
        if self._data is None:
            raise ValueError("self._data must be set, "
                             "using preprocess_data or "
                             "an equivalent method")
        elif self._engine is None:
            raise ValueError("A fitting engine must be set, "
                             "by calling one of: \n"
                             "set_scipy_engine")

        beta = np.max(self._data['y'])*(np.max(self._data['x']) -
                                        np.min(self._data['x']))

        for N in range(1, max_num_features + 1):
            func = self.update_function(func, N)
            self.update_fit_engine(func, params)

            self._engine.do_fit(self._data['x'], self._data['y'],
                                self._data['e'], func)

            params = self.report(func, N, beta)

    def _check_engine_and_data_set_valid(self) -> None:
        """
        A simple check to see if the fit engine
        has already been set and that data is
        available for fitting.
        """
        if self._engine is not None:
            raise RuntimeError("Cannot change the fitting engine.")
        if self._data is None:
            raise ValueError("self._data must be set, "
                             "using preprocess_data or "
                             "an equivalent method")

    def set_scipy_engine(self, guess: ndarray, lower: ndarray,
                         upper: ndarray) -> None:
        """
        Method to set the fit engine to be scipy
        :param guess: the starting guess for the fit
        :param lower: the lower bound for the fit
        :param upper: the upper bound for the fit
        """
        self._check_engine_and_data_set_valid()
        self._engine = ScipyFitEngine(self._data['x'], self._data['y'],
                                      self._data['e'], lower, upper,
                                      guess)

    def update_scipy_fit_engine(self, func: BaseFitFunction, params: ndarray):
        """
        This updates the bounds and guess for scipy fit engine
        :param func: the fitting function
        :param params: the fitting parameters
        """
        lower, upper = func.get_bounds()
        guess = update_guess(list(params), func)
        self._engine.set_guess_and_bounds(guess, lower, upper)
