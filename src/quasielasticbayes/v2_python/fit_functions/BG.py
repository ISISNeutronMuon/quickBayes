from quasielasticbayes.v2.base import BaseFitFunction
from numpy import ndarray
from typing import Dict, List


class LinearBG(BaseFitFunction):
    def __init__(self, prefix=''):
        super().__init__(2, prefix)

    def __call__(self, x: ndarray, m: float, c: float) -> ndarray:
        """
        Implement the Linear BG.
        Need to follow the expected
        form for scipy
        """
        return m*x + c

    def report(self, report_dict: Dict[str, List[float]],
               m: float, c: float) -> Dict[str, List[float]]:
        """
        returns the fit parameters as a dict
        """
        report_dict = self._add_to_report(f"{self._prefix}BG gradient",
                                          m, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}BG constant",
                                          c, report_dict)
        return report_dict
