from quasielasticbayes.v2.base import BaseFitFunction
from numpy import ndarray, pi
import numpy as np
from typing import Dict, List


class Gaussian(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        super().__init__(3, prefix)

    def __call__(self, x: ndarray, amplitude: float, x0: float,
                 sigma: float) -> ndarray:
        """
        Implement a gaussian.
        Need to follow the expected
        form for scipy
        """
        pre_factor = amplitude/(sigma*np.sqrt(2.*pi))
        return pre_factor*np.exp(-pow(x-x0, 2)/(2.*sigma*sigma))

    def report(self, report_dict: Dict[str, List[float]], a: float,
               x0: float, sigma: float) -> Dict[str, List[float]]:
        """
        returns the fit parameters as a dict
        """
        report_dict = self._add_to_report(f"{self._prefix}Amplitude",
                                          a, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Mean",
                                          x0, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Sigma",
                                          sigma, report_dict)
        return report_dict
