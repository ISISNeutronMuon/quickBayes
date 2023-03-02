from quasielasticbayes.v2.functions.composite import CompositeFunction
from quasielasticbayes.v2.functions.exp_decay import ExpDecay
from quasielasticbayes.v2.utils.general import (update_guess,
                                                get_background_function)
from quasielasticbayes.v2.utils.crop_data import crop
from quasielasticbayes.v2.workflow.template import Workflow
from quasielasticbayes.v2.functions.base import BaseFitFunction
from numpy import ndarray
from typing import Dict, List


class MuonExpDecay(Workflow):
    """
    A class for the muon exponential decay workflow
    """
    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        start_x: float, end_x: float) -> None:
        sx, sy, se = crop(x_data, y_data, e_data,
                          start_x, end_x)
        super().preprocess_data(sx, sy, se)

    def _update_function(self, func: BaseFitFunction) -> BaseFitFunction:
        exp_function = ExpDecay()
        func.add_function(exp_function)
        return func

    def update_fit_engine(self, func: BaseFitFunction, params: ndarray):
        lower, upper = func.get_bounds()
        guess = update_guess(list(params), func)
        # assume scipy
        self._engine.set_guess_and_bounds(guess, lower, upper)


def muon_expdecay_main(sample: Dict[str, ndarray],
                       BG_type: str, start_x: float, end_x: float,
                       results: Dict[str, ndarray],
                       results_errors: Dict[str, ndarray],
                       init_params: List[float] = None) -> (Dict[str, ndarray],
                                                            Dict[str, ndarray],
                                                            ndarray,
                                                            List[ndarray],
                                                            List[ndarray]):
    """
    The main function for calculating muon decay rates.
    """
    # construct fitting function
    BG = get_background_function(BG_type)
    func = CompositeFunction()
    func.add_function(BG)
    lower, upper = func.get_bounds()

    # setup workflow
    workflow = MuonExpDecay(results, results_errors)
    workflow.preprocess_data(sample['x'], sample['y'], sample['e'],
                             start_x, end_x)
    params = init_params if init_params is not None else func.get_guess()
    workflow.set_scipy_engine(params, lower, upper)

    # do the calculation
    max_features = 4
    func = workflow.execute(max_features, func, params)
    results, results_errors = workflow.get_parameters_and_errors

    engine = workflow.fit_engine
    fits = []
    errors_fit = []
    x_data = []
    for j in range(max_features):
        x_data, y, e, df, de = engine.get_fit_values(j)
        fits.append(y)
        errors_fit.append(e)

    return results, results_errors, x_data, fits, errors_fit
