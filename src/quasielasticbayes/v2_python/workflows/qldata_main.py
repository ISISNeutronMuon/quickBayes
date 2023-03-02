from quasielasticbayes.v2.functions.qldata_function import QlDataFunction
from quasielasticbayes.v2.utils.spline import spline
from quasielasticbayes.v2.utils.general import (update_guess,
                                                get_background_function)
from quasielasticbayes.v2.workflow.template import Workflow
from quasielasticbayes.v2.functions.base import BaseFitFunction

from numpy import ndarray
import numpy as np
from typing import Dict, List


class QLData(Workflow):
    """
    A class for the quaielastic lorentzians workflow
    """
    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        start_x: float, end_x: float,
                        res: Dict[str, ndarray]) -> (ndarray, ndarray):
        dx = x_data[1] - x_data[0]
        new_x = np.linspace(start_x, end_x, int((end_x - start_x)/dx))

        sy = spline(x_data, y_data, new_x)
        se = spline(x_data, e_data, new_x)
        ry = spline(res['x'], res['y'], new_x)
        super().preprocess_data(new_x, sy, se)

        return new_x, ry

    def _update_function(self, func: BaseFitFunction) -> BaseFitFunction:
        func.add_single_lorentzian()
        return func

    def update_fit_engine(self, func: BaseFitFunction, params: ndarray):
        lower, upper = func.get_bounds()
        guess = update_guess(list(params), func)
        # assume scipy
        self._engine.set_guess_and_bounds(guess, lower, upper)


def ql_data_main(sample: Dict[str, ndarray], res: Dict[str, ndarray],
                 BG_type: str, start_x: float, end_x: float,
                 elastic: bool,
                 results: Dict[str, ndarray],
                 results_errors: Dict[str, ndarray],
                 init_params: List[float] = None) -> (Dict[str, ndarray],
                                                      Dict[str, ndarray],
                                                      ndarray,
                                                      List[ndarray],
                                                      List[ndarray]):
    """
    """
    # setup workflow
    workflow = QLData(results, results_errors)
    new_x, ry = workflow.preprocess_data(sample['x'], sample['y'],
                                         sample['e'],
                                         start_x, end_x, res)

    max_num_peaks = 3

    # setup fit function
    BG = get_background_function(BG_type)
    func = QlDataFunction(BG, elastic, new_x, ry, start_x, end_x)
    lower, upper = func.get_bounds()

    params = init_params if init_params is not None else func.get_guess()
    # just want a guess the same length as lower, it is not used
    workflow.set_scipy_engine(func.get_guess(), lower, upper)

    # do the calculation
    func = workflow.execute(max_num_peaks, func, params)
    results, results_errors = workflow.get_parameters_and_errors

    engine = workflow.fit_engine
    fits = []
    errors_fit = []
    x_data = []
    for j in range(max_num_peaks):
        x_data, y, e, df, de = engine.get_fit_values(j)
        fits.append(y)
        errors_fit.append(e)
    return results, results_errors, x_data, fits, errors_fit
