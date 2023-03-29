from quickBayes.functions.base import BaseFitFunction
from quickBayes.workflow.template import WorkflowTemplate

from quickBayes.log_likelihood import loglikelihood

from numpy import ndarray
import numpy as np
from typing import Dict
from abc import abstractmethod


class ModelSelectionWorkflow(WorkflowTemplate):
    """
    This is a class for the quick bayes model selection workflow.
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
    - set_gofit_engine (gofit)

    Other methods:
    - preprocess_data
    - update_fit_engine
    - update_function (call this one not the overwritten one)
    - report
    - execute
    """

    def __init__(self, results: Dict[str, ndarray],
                 results_errors: Dict[str, ndarray]):
        """
        Set the results and error dicts for reporting
        :param results: dict of parameter values
        :param results_errors: dict of parameter errors
        """
        self._results_dict = results
        self._errors_dict = results_errors
        super().__init__()

    @property
    def get_parameters_and_errors(self) -> (Dict[str, float],
                                            Dict[str, float]):
        """
        Method to get the dict's of the fit parameters and errors.
        :return dict of fit parameters, dict of fit parameter errors
        """
        return self._results_dict, self._errors_dict

    @abstractmethod
    def _update_function(self, func: BaseFitFunction) -> BaseFitFunction:
        """
        This method updates the fitting function when the number of
        features has been incremented.
        It will be unique to the workflow.
        :param func: the fitting function that needs modifying
        :return the modified fitting function
        """
        raise NotImplementedError()

    def update_function(self, func: BaseFitFunction,
                        N: int) -> BaseFitFunction:
        """
        This method updates the fitting function when the number of
        features has been incremented.
        It will be unique to the workflow.
        :param func: the fitting function that needs modifying
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
