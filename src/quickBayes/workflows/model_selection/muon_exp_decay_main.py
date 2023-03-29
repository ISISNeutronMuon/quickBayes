from quickBayes.functions.composite import CompositeFunction
from quickBayes.functions.exp_decay import ExpDecay
from quickBayes.utils.general import get_background_function
from quickBayes.utils.crop_data import crop
from quickBayes.workflow.model_template import ModelSelectionWorkflow
from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray
from typing import Dict, List


class MuonExpDecay(ModelSelectionWorkflow):
    """
    A class for the muon exponential decay workflow
    """
    def preprocess_data(self, x_data: ndarray,
                        y_data: ndarray, e_data: ndarray,
                        start_x: float, end_x: float) -> None:
        """
        The preprocessing needed for the data.
        This crops and stores the data.
        :param x_data: the x data to fit to
        :param y_data: the y data to fit to
        :param e_data: the errors for the y data
        :param start_x: the start x value
        :param end_x: the end x value
        """
        sx, sy, se = crop(x_data, y_data, e_data,
                          start_x, end_x)
        super().preprocess_data(sx, sy, se)

    @staticmethod
    def _update_function(func: BaseFitFunction) -> BaseFitFunction:
        """
        This method adds a exponential decay to the fitting
        function.
        :param func: the fitting function that needs modifying
        :return the modified fitting function
        """

        exp_function = ExpDecay()
        func.add_function(exp_function)
        return func


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
    Uses the muon exp decay workflow
    :param sample: dict containing the sample x, y and e data (keys = x, y, e)
    :param BG_type: the type of BG ("none", "flat", "linear")
    :param start_x: the start x for the calculation
    :param end_x: the end x for the calculation
    :param results: dict of results
    :param results_errors: dict of errors for results
    :param init_params: initial values, if None a guess will be made
    :result dict of the fit parameters, their errors, the x range used, list of
    fit values and their errors.
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
