from quasielasticbayes.v2.workflow.template import Workflow
from quasielasticbayes.v2.log_likelihood import loglikelihood

from numpy import ndarray
import numpy as np
from typing import Dict
from abc import abstractmethod


class Axis(object):
    def __init__(self, start, end, N, label='x'):
        self._len = N
        self._label = label
        self._vals = np.linspace(start, end, N)

    @property
    def label(self):
        return self._label

    @property
    def len(self):
        return self._len

    @property
    def values(self):
        return self._vals


class GridSearchTemplate(Workflow):
    def __init__(self, results: Dict[str, ndarray],
                 results_errors: Dict[str, ndarray]):
        """
        Set the results and error dicts for reporting
        :param results: dict of parameter values
        :param results_errors: dict of parameter errors
        """
        super().__init__(results, results_errors)
        self._x_axis = None
        self._y_axis = None
        self._grid = None

    def set_x_axis(self, start, end, N):
        self._x_axis = Axis(start, end, N)

    @property
    def get_x_axis(self):
        return self._x_axis

    def set_y_axis(self, start, end, N):
        self._y_axis = Axis(start, end, N)

    @property
    def get_y_axis(self):
        return self._y_axis

    def generate_mesh(self):
        if self._x_axis is None or self._y_axis is None:
            raise RuntimeError("The x and/or y axis has"
                               " not been set. Please use "
                               "set_x_axis and set_y_axis.")

        X, Y = np.meshgrid(self._x_axis.values,
                           self._y_axis.values)
        self._grid = self.empty_mesh(X, Y)
        return X, Y

    @staticmethod
    def empty_mesh(X, Y):
        z = X + Y
        return np.zeros(z.shape)

    @property
    def get_grid(self):
        return self._grid

    def get_slices(self):
        indicies = np.where(self._grid == self._grid.max())
        x_slice = self._grid[indicies[0], :]
        y_slice = self._grid[:, indicies[1]]
        return x_slice, y_slice

    @abstractmethod
    def set_x_value(self, func, value):
        raise NotImplementedError()

    @abstractmethod
    def set_y_value(self, func, value):
        raise NotImplementedError()

    def get_z_value(self, N_x, N_peaks, scale):
        chi2 = self._engine.get_chi_squared()
        covar = self._engine.get_covariance_matrix()
        return loglikelihood(N_x, chi2, covar,
                             N_peaks, scale)

    def normalise_grid(self):
        norm = np.min(self._grid)/self._grid
        norm -= np.min(self._grid)
        norm /= np.max(self._grid)
        self._grid = norm

    @abstractmethod
    def N(self):
        raise NotImplementedError()

    def execute(self, func, results):
        x_data = self._data['x']
        y_data = self._data['y']
        e_data = self._data['e']
        scale = np.max(y_data)*(np.max(x_data) - np.min(x_data))

        self.generate_mesh()

        N = len(self.get_x_axis.values)*len(self.get_y_axis.values)
        counter = 0

        for i, xx in enumerate(self.get_x_axis.values):
            func = self.set_x_value(func, xx)
            params = func.get_guess()

            for j, yy in enumerate(self.get_y_axis.values):
                func = self.set_y_value(func, yy)
                self._engine.do_fit(x_data, y_data, e_data, func)
                params, _ = self._engine.get_fit_parameters()

                results = func.report(results, params)
                self._grid[j][i] = self.get_z_value(len(x_data),
                                                    self.N,
                                                    scale)
                counter += 1
                print(f'\rPercentage complete: {100*counter/N:2f}')
        self.normalise_grid()
