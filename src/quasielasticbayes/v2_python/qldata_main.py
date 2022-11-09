from quasielasticbayes.v2.functions.qldata_function import QlDataFunction
from quasielasticbayes.v2.fitting.scipy_fit import scipy_curve_fit
from quasielasticbayes.v2.utils.spline import spline
from quasielasticbayes.v2.utils.general import (update_guess,
                                                get_background_function)
from quasielasticbayes.v2.log_likelihood import loglikelihood

from numpy import ndarray
import numpy as np
from typing import Dict, List


def ql_data_main(sample: Dict[str, ndarray], res: Dict[str, ndarray],
                 BG_type: str, start_x: float, end_x: float,
                 elastic: bool) -> (Dict[str, ndarray], List[float]):
    """
    The main function for calculating Qldata.
    Steps are:
    1. Interpolate and crop the data to have same bins
    2. Fit BG + delta (if elastic=True) + lorentzian
    3. Calculate loglikelihood for fit
    4. Fit BG + delta (if elastic=True) + 2 lorentzians
    5. Calculate loglikelihood for fit
    6. Fit BG + delta (if elastic=True) + 3 lorentzians
    7. Calculate loglikelihood for fit
    8. Report results
    :param sample: dict containing the sample x, y and e data (keys = x, y, e)
    :param res: dict containg the resolution x, y data (keys = x, y)
    :param BG_type: the type of BG ("none", "flat", "linear")
    :param start_x: the start x for the calculation
    :param end_x: the end x for the calculation
    :param elastic: if to include the elastic peak
    :result dict of the fit parameters and an array of the loglikelihoods
    """
    # step 0
    results = {}
    BG = get_background_function(BG_type)
    max_num_peaks = 3

    # step 1
    x_data = sample['x']
    dx = x_data[1] - x_data[0]
    new_x = np.linspace(start_x, end_x, int((end_x - start_x)/dx))

    sy = spline(sample['x'], sample['y'], new_x)
    se = spline(sample['x'], sample['e'], new_x)
    ry = spline(res['x'], res['y'], new_x)

    func = QlDataFunction(BG, elastic, new_x, ry, start_x, end_x)
    guess = func.get_guess()
    # loop doing steps 2 to 8
    params = guess
    for N in range(1, max_num_peaks+1):
        func.add_single_lorentzian()

        lower, upper = func.get_bounds()
        params = update_guess(params, func)
        (chi2, hess_det, params, fit) = scipy_curve_fit(new_x, sy, se,
                                                        func, params,
                                                        lower, upper)
        params = list(params)

        results = func.report(results, *params)

        prob_name = 'N{N}:loglikelihood'
        if prob_name in results:
            results[prob_name].append(loglikelihood(len(sy), chi2,
                                                    hess_det,
                                                    func.N_peaks))
        else:
            results[f'N{N}:loglikelihood'] = [loglikelihood(len(sy), chi2,
                                                            hess_det,
                                                            func.N_peaks)]

    return results, new_x
