"""
Implements a controller for the Bumps fitting software.
"""

from bumps.fitters import fit as bumpsFit
from bumps.names import Curve, FitProblem, PoissonCurve

from quickBayes.fitting.fit_engine import FitEngine
import numpy as np

import scipy
import numdifftools as ndt
from quickBayes.fitting.fit_utils import (chi_squared,
                                          param_errors,
                                          derivative,
                                          fit_errors,
                                          var, res)

def wrapper(exec_dict, names):
    _names = ','.join(names)
    wrap = f'def fitFunction(x, {_names}):\n'
    #wrap += f'    return func([{_names}], x=x)'
    wrap += f'    return func(x, [{_names}])'
    exec(wrap, exec_dict)
    return wrap

class Bumps(FitEngine):

    def __init__(self, func, x, y, e):
        super().__init__('bumps', x, y, e)
        self._func = func
        self._result = None
        self._fit_problem = None
        self._minimizer = 'scipy.leastsq'

    def set_cost_function(self, func):
        self.cost_function = func

    def cost_function(self, x, p):
        tmp = (self._func(x, *p) - self._y_data)/self._e_data
    
        return np.ravel(tmp)#np.sum(tmp)

    #    #return res(self._func, x, self._y_data, self._e_data, p)

    def set_function(self, function, full_names=[]):
        # pylint: disable=exec-used,protected-access
        """
        Setup problem ready to run with Bumps.

        Creates a FitProblem for calling in the fit() function of Bumps
        """
        self._func = function
        # create some dummy param names
        #if params == []:
        
        guess = function.get_guess()
        #else:
        #guess = params
        param_names = full_names if full_names != [] else [ f'a{j}' for j in range(len(guess))]
        # Bumps fails with the *args notation
        param_name_str = ', '.join(param_names)
        param_dict = dict(zip(param_names, guess))

        

        # Create a Function Wrapper for the problem function. The type of the
        # Function Wrapper is acceptable by Bumps.
        # Send in the residual as the model, with zero
        # y data.  This allows all our supported nlls
        # cost fucntions to be used.
        #exec_dict = {'func': self.cost_function}
        exec_dict = {'func': self.cost_function}
        #wrapper = self._wrapper(exec_dict, param_name_str)
        wrap = wrapper(exec_dict, param_names)
        model = exec_dict['fitFunction']
        func_wrapper = Curve(fn=model,
                                 x=self._x_data,
                                 y=np.zeros(len(self._y_data)),
                                 **param_dict)

        # Set a range for each parameter
        lower, upper = self._func.get_bounds()
        for ind, name in enumerate(param_names):
            #if self.value_ranges is not None:
            #    min_val = l[ind]
            #    max_val = u[ind]
            #else:
            min_val = lower[ind]
            max_val = upper[ind]
            func_wrapper.__dict__[name].range(min_val, max_val)

        # Create a Problem Wrapper. The type of the Problem Wrapper is
        # acceptable by Bumps fitting.
        self._func_wrapper = func_wrapper
        self._fit_problem = FitProblem(func_wrapper)
        self._minimizer = "scipy.leastsq"


        # Determine the order of the parameters in `self.fit_problem` as this
        # could differ from the ordering of parameters in `self._param_names`
        param_order = []
        for i in range(len(param_names)):
            param_order.append(str(self._fit_problem._parameters[i]))
        self.fit_order = param_order

    
    """
    def do_fit(self, x_data, y_data, e_data
               , func):
        params = self._do_fit(x_data, y_data, e_data, func)

        def ff(p):
            return self._func(self._x_data, *p)
        def rr(p):
            return var(func, self._x_data, self._y_data, p)

        tmp = ndt.Gradient(ff)(params)
        df_by_dp = tmp.T # derivative(x_data, params, func)
        H = ndt.Hessian(rr, step=1e-4)
        hess = H(params)
        

        self._covars.append(scipy.linalg.inv(hess)*2.) #calculate_covar(x_data, y_data, e_data, func, df_by_dp, params)
        self.add_params(params)
        self.add_fit(self._x_data, self._func, df_by_dp, params)
        self._chi2.append(chi_squared(self._x_data, self._y_data, self._e_data,
                                      self._fit, params))
        self._fit = None
    """
    def _do_fit(self, x,y,e, func):
        result = self.fit()
        return result.x

    def fit(self):
        """
        Run problem with Bumps.
        """
        return bumpsFit(self._fit_problem,
                        method=self._minimizer)
