from quickBayes.fitting.scipy_engine import ScipyFitEngine
from quickBayes.fitting.gofit_engine import GoFitEngine
from quickBayes.functions.base import BaseFitFunction

from quickBayes.utils.general import update_guess

from numpy import ndarray
from abc import abstractmethod


class WorkflowTemplate(object):
    """
    This is a class for the quick bayes workflow.
    Each method can be overwritten to provide
    unique functionality for the specific
    use case.

    The properties are:
    - fit_engine

    To add a fit engine:
    - set_scipy_engine (scipy curve fit)
    - set_gofit_engine (gofit)

    Other methods:
    - preprocess_data
    - update_fit_engine
    - update_function (call this one not the overwritten one)
    - execute
    """
    def __init__(self):
        """
        Set the results and error dicts for reporting
        """
        self._engine = None
        self._data = None
        self._raw = None

    @property
    def get_raw(self):
        """
        Returns the original data
        (i.e. the data fits are interpolated to match)
        :return a dict of the raw data
        """
        return self._raw

    @property
    def fit_engine(self):
        """
        Simple method for getting the fit engine
        :return the fit engine used
        """
        return self._engine

    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        *args: float) -> None:
        """
        The preprocessing needed for the data.
        This simple case just assigns the data values.
        This can be overwritten in a derived class to
        include extra processing (e.g. cropping or
        interpolating the data).
        The data and raw data are assumed to be the same here
        :param x_data: the x data to fit to
        :param y_data: the y data to fit to
        :param e_data: the errors for the y data
        :param *args: additional arguments that might be needed
        """
        self._data = {'x': x_data, 'y': y_data, 'e': e_data}
        self._raw = {'x': x_data, 'y': y_data, 'e': e_data}

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
        elif self._engine.name == 'gofit':
            self.update_gofit_engine(func)
        else:
            raise RuntimeError("The fit engine is "
                               f"{self._engine.name} "
                               "please use the appropriate "
                               "update method")
        return

    @abstractmethod
    def execute(self, *args) -> None:
        """
        The main part of the analysis.
        It increments the number of features in the fitting function,
        does a fit and then records the results.
        :param args: the arguments
        """
        raise NotImplementedError()

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
        self._engine = ScipyFitEngine(self._raw['x'], self._raw['y'],
                                      self._raw['e'], lower, upper,
                                      guess)

    def _get_bounds(self, func: BaseFitFunction) -> (ndarray, ndarray):
        """
        Get the bounds for the fit engine
        :param func: the fit function
        :returns the lower and upper bounds
        """
        return func.get_bounds()

    def update_scipy_fit_engine(self, func: BaseFitFunction, params: ndarray):
        """
        This updates the bounds and guess for scipy fit engine
        :param func: the fitting function
        :param params: the fitting parameters
        """
        lower, upper = self._get_bounds(func)

        guess = update_guess(list(params), func)
        self._engine.set_guess_and_bounds(guess, lower, upper)

    def set_gofit_engine(self, samples: int, lower: ndarray,
                         upper: ndarray) -> None:
        """
        Method to set the fit engine to be gofit
        :param samples: the number of samples to use
        :param lower: the lower bound for the fit
        :param upper: the upper bound for the fit
        """
        self._check_engine_and_data_set_valid()
        self._engine = GoFitEngine(self._raw['x'], self._raw['y'],
                                   self._raw['e'], lower, upper, samples)

    def update_gofit_engine(self, func: BaseFitFunction):
        """
        This updates the bounds for gofit engine
        :param func: the fitting function
        """
        lower, upper = func.get_bounds()
        self._engine.set_bounds_and_N_params(lower, upper)
