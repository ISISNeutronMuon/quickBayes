from quasielasticbayes.v2.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List


class Delta(BaseFitFunction):
    """
    Strictly this is not a true delta function.
    Instead it is a top hat function, which
    in the limit of binwidth-> 0 is a delta
    """
    def __init__(self, prefix: str = ''):
        super().__init__(2, prefix)

    def __call__(self, x: ndarray, amplitude: float, x0: float) -> ndarray:
        """
        Implement the delta/top hat.
        Need to follow the expected
        form for scipy
        """
        data = np.zeros(len(x))
        index = np.searchsorted(x, x0)-1

        # integral should normalise to 1*amplitude
        # so need to divide by bin width
        dx = 0.0
        if index == len(data)-1:
            dx = x[index] - x[index-1]
        else:
            dx = x[index+1] - x[index]

        data[index] = amplitude/dx
        return data

    def report(self, report_dict: Dict[str, List[float]],
               amplitude: float, x0: float) -> Dict[str, List[float]]:
        """
        returns the fit parameters as a dict
        """
        report_dict = self._add_to_report(f"{self._prefix}Amplitude",
                                          amplitude, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Centre",
                                          x0, report_dict)
        return report_dict
