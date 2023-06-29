from quickBayes.workflow.template import WorkflowTemplate
from quickBayes.log_likelihood import loglikelihood
from quickBayes.functions.base import BaseFitFunction

from numpy import ndarray
import numpy as np
from abc import abstractmethod


class Axis(object):
    """
    A simple object for holding the information
    about an axis
    """
    def __init__(self, start, end, N, label='x'):
        """
        Create an axis
        :param start: the first value for the axis
        :param end: the last value for the axis
        :param N: the number of values along the axis
        :param label: the name of the axis
        """
        self._len = N
        self._label = label
        self._vals = np.linspace(start, end, N)

    @property
    def label(self) -> str:
        """
        Get the label for the axis
        :return the label
        """
        return self._label

    @property
    def len(self):
        """
        Get the length of the axis
        :return the length of the axis
        """
        return self._len

    @property
    def values(self):
        """
        Get the values used on the axis
        :return the values from the axis
        """
        return self._vals


class GridSearchTemplate(WorkflowTemplate):
    """
    A workflow for a grid search.
    A grid search will do a series of fits for a range of
    x and y values, that correspond to fixed parameters in
    the fitting function.

    The inherited class must include:
    - _set_x_value
    - _set_y_value
    - N

    The properties are:
    - fit_engine
    - get_grid
    - get_x_axis
    - get_y_axis

    To add a fit engine:
    - set_scipy_engine (scipy curve fit, recommended)
    - set_gofit_engine (gofit)

    Other methods:
    - preprocess_data
    - update_fit_engine
    - update_function (call this one not the overwritten one)
    - execute
    - set_x_axis
    - set_y_axis
    - N
    """
    def __init__(self):
        """
        Set the results and error dicts for reporting
        """
        super().__init__()
        self._x_axis = None
        self._y_axis = None
        self._grid = None

    def set_x_axis(self, start: float, end: float,
                   N: int, label: str) -> None:
        """
        Sets the x axis for the workflow
        :param start: the first value on the x axis
        :param end: the last value on the x axis
        :param N: the number of values on the x axis
        :param label: the x axis label
        """
        self._x_axis = Axis(start, end, N, label)

    @property
    def get_x_axis(self) -> Axis:
        """
        Get the x axis
        :return the x axis object
        """
        return self._x_axis

    def set_y_axis(self, start: float, end: float,
                   N: int, label: str) -> None:
        """
        Sets the y axis for the workflow
        :param start: the first value on the y axis
        :param end: the last value on the y axis
        :param N: the number of values on the y axis
        :param label: the y axis label
        """
        self._y_axis = Axis(start, end, N, label)

    @property
    def get_y_axis(self) -> Axis:
        """
        Get the y axis
        :return the y axis object
        """
        return self._y_axis

    def _generate_grid(self) -> (ndarray, ndarray):
        """
        Creates a grid from the axis.
        Need to have set the x and y axis first.
        :return the X and Y values from np.meshgrid
        """
        if self._x_axis is None or self._y_axis is None:
            raise ValueError("The x and/or y axis has"
                             " not been set. Please use "
                             "set_x_axis and set_y_axis.")

        X, Y = np.meshgrid(self._x_axis.values,
                           self._y_axis.values)
        self._grid = self._empty_mesh(X)
        return X, Y

    @staticmethod
    def _empty_mesh(X: ndarray) -> None:
        """
        Creates a grid of the correct size with zeros
        :param X: one of the meshgrid outputs
        """
        return np.zeros(X.shape)

    @property
    def get_grid(self) -> ndarray:
        """
        Get the grid values
        :return the grid values
        """
        return self._grid

    def get_slices(self) -> (ndarray, ndarray):
        """
        Gets slices along the x and y axis, such that
        the peak grid value is in both slices.
        :return the x and y data slices
        """
        indices = np.where(self._grid == self._grid.max())
        x_slice = self._grid[indices[0], :][0]
        y_slice = [vec[0] for vec in self._grid[:, indices[1]]]
        return x_slice, np.array(y_slice)

    @abstractmethod
    def _set_x_value(self, func: BaseFitFunction,
                     value: float) -> BaseFitFunction:
        """
        Sets the x value (fixed fit parameter)
        :param func: the function that is being updated
        :param value: the fix value
        """
        raise NotImplementedError()

    @abstractmethod
    def _set_y_value(self, func: BaseFitFunction,
                     value: float) -> BaseFitFunction:
        """
        Sets the y value (fixed fit parameter)
        :param func: the function that is being updated
        :param value: the fix value
        """
        raise NotImplementedError()

    def _get_z_value(self, N_p: int, N_f: int, scale: float) -> float:
        """
        Gets the loglikelihood value
        :param N_p: the number of data points being fitted to
        :param N_f: the number of features being fitted to
        :param scale: the beta scale factor for loglikelihood
        :return the loglikelihood
        """
        chi2 = self._engine.get_chi_squared()
        covar = self._engine.get_covariance_matrix()
        return loglikelihood(N_p, chi2, covar,
                             N_f, scale)

    def _normalise_grid(self) -> None:
        """
        Rescales the grid values so that the max is 1
        and the min is 0. This should make features clearer
        """
        self._grid = np.min(self._grid)/self._grid
        self._grid -= np.min(self._grid)
        self._grid /= np.max(self._grid)

    @abstractmethod
    def N(self, func: BaseFitFunction) -> int:
        """
        Gets the number of features in the fit function
        :param func: the fitting function
        :return the number of features
        """
        raise NotImplementedError()

    def _get_bounds(self, func: BaseFitFunction) -> (ndarray, ndarray):
        """
        Get the bounds for the fit engine
        For a grid search we want the original
        bounds
        :param func: the fit function
        :returns the lower and upper bounds
        """
        return self._engine._lower, self._engine._upper

    def execute(self, func: BaseFitFunction) -> (ndarray, ndarray):
        """
        Does the grid search. Needs the x and y axis to be set.
        Also needs a fitting engine to be set.
        :param func: the fitting function
        :return the X and Y values for the grid
        """
        if self._engine is None:
            raise ValueError("please set a fit engine")
        x_data = self._data['x']
        y_data = self._data['y']
        e_data = self._data['e']
        scale = np.max(y_data)*(np.max(x_data) - np.min(x_data))

        X, Y = self._generate_grid()
        for i, xx in enumerate(self.get_x_axis.values):
            func = self._set_x_value(func, xx)
            params = func.get_guess()

            for j, yy in enumerate(self.get_y_axis.values):
                func = self._set_y_value(func, yy)
                self._engine.do_fit(x_data, y_data, e_data, func)
                params, _ = self._engine.get_fit_parameters()

                num = self.N(func)
                self._grid[j][i] = self._get_z_value(len(x_data),
                                                     num,
                                                     scale)
                self.update_fit_engine(func, params)
        self._normalise_grid()
        return X, Y
