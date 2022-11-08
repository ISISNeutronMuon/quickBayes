from quasielasticbayes.v2.functions.base import BaseFitFunction
from numpy import ndarray
from typing import Dict, List


class LinearBG(BaseFitFunction):
    def __init__(self, prefix=''):
        """
        :param prefix: prefix for function parameters in report
        """
        super().__init__(2, prefix)

    def __call__(self, x: ndarray, m: float, c: float) -> ndarray:
        """
        Implement the Linear BG.
        Need to follow the expected
        form for scipy
        :param x: x values
        :param m: gradient
        :param c: constant
        :return linear background y values
        """
        return m*x + c

    def report(self, report_dict: Dict[str, List[float]],
               m: float, c: float) -> Dict[str, List[float]]:
        """
        reporting method
        :param report_dict: dict of parameters
        :param m: gradient
        :param c: constant
        :return dict of parameters, including BG
        """
        report_dict = self._add_to_report(f"{self._prefix}BG gradient",
                                          m, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}BG constant",
                                          c, report_dict)
        return report_dict

    def get_guess(self) -> List[float]:
        """
        Get the starting guess for a fit function
        :return the initial guess
        """
        return [0., 0.]

    def get_bounds(self) -> (List[float], List[float]):
        """
        Get the fitting bounds
        :return lists for lower and upper bounds
        """
        return [-1., -1.], [1., 1.]
