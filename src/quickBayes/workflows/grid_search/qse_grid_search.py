from quickBayes.workflow.grid_template import GridSearchTemplate
from quickBayes.functions.qse_fixed import QSEFixFunction
from quickBayes.utils.spline import spline
from numpy import ndarray
from typing import Dict
import numpy as np


class QSEGridSearch(GridSearchTemplate):
    """
    A workflow for doing a grid search of quasielastic
    data to fit a stretch exponential for different
    (fixed) beta and FWHM values.

    The properties are:
    - fit_engine
    - get_grid
    - get_x_axis
    - get_y_axis
    - N

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
    """
    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        start_x: float, end_x: float,
                        res: Dict[str, ndarray]) -> (ndarray, ndarray):
        """
        The preprocessing needed for the data.
        It splines the sample and resolution data
        to the same uniform grid.
        This is designed for a fixed stretched exp.
        :param x_data: the sample x data to fit to
        :param y_data: the sample y data to fit to
        :param e_data: the sample errors for the y data
        :param start_x: the start x value
        :param end_x: the end x value
        :param res: a dict of the resolution data (keys =x, y, e)
        :return the new x range and the new resolution y values
        """
        dx = x_data[1] - x_data[0]
        new_x = np.linspace(start_x, end_x, int((end_x - start_x)/dx))

        sy = spline(x_data, y_data, new_x)
        se = spline(x_data, e_data, new_x)
        ry = spline(res['x'], res['y'], new_x)
        super().preprocess_data(new_x, sy, se)

        return new_x, ry

    @staticmethod
    def _set_x_value(func: QSEFixFunction,
                     value: float) -> QSEFixFunction:
        """
        Sets the beta value for the fit
        function (x axis)
        :param func: the stretch exp with fixes function
        :param value: the value to fix beta to
        :return the updated fit function
        """
        func.set_beta(value)
        return func

    @staticmethod
    def _set_y_value(func: QSEFixFunction,
                     value: float) -> QSEFixFunction:
        """
        Sets the FWHM (tau) value for the fit
        function (y axis)
        :param func: the stretch exp with fixes function
        :param value: the value to fix FWHM (tau) to
        :return the updated fit function
        """
        func.set_FWHM(value)
        return func

    @staticmethod
    def N(func: QSEFixFunction) -> int:
        """
        Get the number of features from fit
        :param func: the fitting function
        :return the number of features
        """
        return func.N_peaks
