from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List


class NoBG(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        :param prefix: prefix for function parameters in report
        """
        super().__init__(0, prefix, [], [], [])

    def __call__(self, x: ndarray) -> ndarray:
        """
        Implement mo BG.
        Need to follow the expected
        form for scipy
        :param x: x values
        :return no background
        """
        return np.zeros(len(x))

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return []

    def report(self, report_dict: Dict[str,
                                       List[float]]) -> Dict[str, List[float]]:
        """
        reporting method
        :param report_dict: dict of parameters
        :return dict of parameters
        """
        return report_dict


class FlatBG(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        :param prefix: prefix for function parameters in report
        """
        super().__init__(1, prefix, [0.], [-1.], [1.])

    @property
    def constant(self) -> str:
        return str(f'{self._prefix}BG constant')

    def __call__(self, x: ndarray, c: float) -> ndarray:
        """
        Implement the flat BG.
        Need to follow the expected
        form for scipy
        :param x: x values
        :param c: constant
        :return linear background y values
        """
        return c*np.ones(len(x))

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.constant, index)]

    def report(self, report_dict: Dict[str, List[float]],
               c: float) -> Dict[str, List[float]]:
        """
        reporting method
        :param report_dict: dict of parameters
        :param c: constant
        :return dict of parameters, including BG
        """
        report_dict = self._add_to_report(self.constant,
                                          c, report_dict)
        return report_dict


class LinearBG(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        :param prefix: prefix for function parameters in report
        """
        super().__init__(2, prefix, [0., 0.], [-1., -1.], [1., 1.])

    @property
    def constant(self) -> str:
        return str(f'{self._prefix}BG constant')

    @property
    def grad(self) -> str:
        return str(f'{self._prefix}BG gradient')

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

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.grad, index),
                self._read_report(report_dict, self.constant, index)]

    def report(self, report_dict: Dict[str, List[float]],
               m: float, c: float) -> Dict[str, List[float]]:
        """
        reporting method
        :param report_dict: dict of parameters
        :param m: gradient
        :param c: constant
        :return dict of parameters, including BG
        """
        report_dict = self._add_to_report(self.grad,
                                          m, report_dict)
        report_dict = self._add_to_report(self.constant,
                                          c, report_dict)
        return report_dict
