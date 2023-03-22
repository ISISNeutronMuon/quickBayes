from quasielasticbayes.v2.workflow.grid_template import GridSearchTemplate
# from quasielasticbayes.v2.functions.qse_fixed import QSEFunction
from quasielasticbayes.v2.utils.spline import spline
from numpy import ndarray
import numpy as np
from typing import Dict


class QSEGridSearch(GridSearchTemplate):

    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        start_x: float, end_x: float,
                        res: Dict[str, ndarray]) -> (ndarray, ndarray):
        """
        The preprocessing needed for the data.
        It splines the sample and resolution data
        to the same uniform grid.
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
    def set_x_value(func, value):
        func.set_beta(value)
        return func

    @staticmethod
    def set_y_value(func, value):
        func.set_FWHM(value)
        return func

    @staticmethod
    def N():
        return 1
