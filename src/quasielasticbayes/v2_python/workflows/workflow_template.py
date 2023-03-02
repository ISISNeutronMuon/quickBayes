from quasielasticbayes.v2.fitting.scipy_engine import ScipyFitEngine
from quasielasticbayes.v2.functions.base import BaseFitFunction

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
    """
    def __init__(self, results: Dict[str, ndarray],
                 results_errors: Dict[str, ndarray]):
        self._engine = None
        self._results_dict = results
        self._errors_dict = results_errors
        self._data = None

    @property
    def fit_engine(self):
        return self._engine

    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        *args: float) -> None:
        self._data = {'x': x_data, 'y': y_data, 'e': e_data}

    @abstractmethod
    def _update_function(self, func: BaseFitFunction) -> BaseFitFunction:
        raise NotImplementedError()

    def update_function(self, func: BaseFitFunction,
                        N: int) -> BaseFitFunction:
        function = self._update_function(func)
        function.update_prefix(f'N{N}:')
        return function

    def update_fit_engine(self, func: BaseFitFunction, params: ndarray):
        return

    def _check_engine_and_data_set_valid(self):
        if self._engine is not None:
            raise RuntimeError("Cannot change the fitting engine.")
        if self._data is None:
            raise ValueError("self._data must be set, "
                             "using preprocess_data or "
                             "an equivalent method")

    def set_scipy_engine(self, guess: ndarray, lower: ndarray,
                         upper: ndarray) -> None:
        self._check_engine_and_data_set_valid()
        self._engine = ScipyFitEngine(self._data['x'], self._data['y'],
                                      self._data['e'], lower, upper,
                                      guess)

    def report(self, func: BaseFitFunction, N: int, beta: float) -> None:
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

    @property
    def get_parameters_and_errors(self) -> (Dict[str, float],
                                            Dict[str, float]):
        return self._results_dict, self._errors_dict
