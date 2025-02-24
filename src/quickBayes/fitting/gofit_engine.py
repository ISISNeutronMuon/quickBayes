from numpy import ndarray
from quickBayes.fitting.fit_engine import FitEngine


try:
    from quickBayes.fitting.gofit_base import _GoFitEngine

    class GoFitEngine(_GoFitEngine):

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
            super().__init__(x_data, y_data, e_data,
                             lower, upper, samples, max_iterations)

except ImportError:

    class GoFitEngine(FitEngine):

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
            raise RuntimeError("gofit is not installed. Please "
                               "install gofit to use this functionality")
