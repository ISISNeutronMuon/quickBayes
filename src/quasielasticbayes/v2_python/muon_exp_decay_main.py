from quasielasticbayes.v2.functions.composite import CompositeFunction
from quasielasticbayes.v2.functions.exp_decay import ExpDecay
from quasielasticbayes.v2.fitting.scipy_fit import scipy_curve_fit
from quasielasticbayes.v2.fitting.fit_utils import (log10_hessian_det,
                                                    chi_squared,
                                                    param_errors,
                                                    derivative,
                                                    fit_errors)
from quasielasticbayes.v2.utils.general import (update_guess,
                                                get_background_function)
from quasielasticbayes.v2.log_likelihood import loglikelihood
from numpy import ndarray
import numpy as np
from typing import Dict, List


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
    # step 0
    BG = get_background_function(BG_type)
    max_num = 4
    # step 1
    x_data = sample['x']
    sy = sample['y']
    se = sample['e']

    func = CompositeFunction()
    func.add_function(BG)
    fits = []
    errors_fit = []
    params = init_params if init_params is not None else []
    beta = np.max(sy)*(np.max(x_data) - np.min(x_data))
    # loop doing steps 2 to 8
    for N in range(1, max_num+1):
        exp_function = ExpDecay()
        func.add_function(exp_function)
        func.update_prefix(f'N{N}:')

        lower, upper = func.get_bounds()
        params = update_guess(params, func)
        (params, covar, fit) = scipy_curve_fit(x_data, sy, se,
                                               func, params,
                                               lower, upper)
        fits.append(fit)

        chi2 = chi_squared(x_data, sy, se, fit, params)
        hess_det = log10_hessian_det(covar)

        errors_p = param_errors(covar)
        df_by_dp = derivative(x_data, params, func)
        tmp = fit_errors(x_data, params, fit, covar, df_by_dp)
        errors_fit.append(tmp)

        params = list(params)
        results = func.report(results, *params)
        results_errors = func.report_errors(results_errors, errors_p, params)

        prob_name = f'N{N}:loglikelihood'
        if prob_name in results:
            results[prob_name].append(loglikelihood(len(sy), chi2,
                                                    hess_det,
                                                    N, beta))
        else:
            results[prob_name] = [loglikelihood(len(sy), chi2,
                                                hess_det,
                                                N, beta)]
    return results, results_errors, x_data, fits, errors_fit
