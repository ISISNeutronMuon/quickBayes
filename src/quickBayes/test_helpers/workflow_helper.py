from quickBayes.functions.BG import LinearBG
from quickBayes.functions.composite import CompositeFunction
import numpy as np


"""
Helper functions and classes
for testing workflows
"""

"""
data generation
"""


def gen_model_selection_data():
    x = np.linspace(1, 4, 100)
    np.random.seed(1)
    noise_stdev = 0.1
    y = np.random.normal(0.5 +
                         1.2*np.exp(-2*x), noise_stdev)
    e = 0.5*noise_stdev*np.ones(x.shape)
    return x, y, e


def gen_grid_search_data():
    x = np.linspace(1, 4, 100)
    np.random.seed(1)
    noise_stdev = 0.1
    y = np.random.normal(0.5 + 0.2*x +
                         1.2*np.exp(-2*x), noise_stdev)
    e = 0.5*noise_stdev*np.ones(x.shape)
    return x, y, e


"""
Fitting function, need some definitions
of fitting functions with fixes for
testing grid search.
These are minimal fitting functions
so that the tests can pass. Some
methods/functionality is missing.
"""


class FixedBG(LinearBG):
    def __init__(self, prefix=''):
        super().__init__(prefix)
        self._N_params = 0
        self._guess = []
        self._lower = []
        self._upper = []
        self._c = 0
        self._m = 0

    def set_c(self, val):
        self._c = val

    def set_m(self, val):
        self._m = val

    def report(self, result):
        return super().report(result, self._m, self._c)

    def __call__(self, x):
        return super().__call__(x, self._m, self._c)


class FixedComposite(CompositeFunction):
    def set_c(self, val):
        # assume first entry is always fixed func
        self._funcs[0].set_c(val)

    def set_m(self, val):
        # assume first entry is always fixed func
        self._funcs[0].set_m(val)
