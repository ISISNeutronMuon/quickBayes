from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List


class ExpDecay(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Create an exponential decay function
        :param prefix: prefix for parameter reporting
        """
        super().__init__(2, prefix, [1., 0.1], [0., 0.001], [1., 20.])

    @property
    def amplitude(self) -> str:
        """
        :retrun string for amplitude name
        """
        return str(f"{self._prefix}Amplitude")

    @property
    def decay_rate(self) -> str:
        """
        :return string for name of lambda
        """
        return str(f"{self._prefix}lambda")

    def __call__(self, x: ndarray, amplitude: float,
                 decay_rate: float) -> ndarray:
        """
        Implement an exponential decay.
        Need to follow the expected
        form for scipy
        :param x: x values for the function evaluation
        :param amplitude: amplitude of decay
        :param decay_rate: the lambda value (decay rate)
        :return y values for the function
        """
        return amplitude*np.exp(-decay_rate*x)

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.decay_rate, index)]

    def report(self, report_dict: Dict[str, List[float]], a: float,
               decay_rate: float) -> Dict[str, List[float]]:
        """
        Report the results
        :param report_dict: dict of results
        :param a: amplitude of exp decay
        :param decay_rate: the lambda value
        returns the updated results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          a, report_dict)
        report_dict = self._add_to_report(self.decay_rate,
                                          decay_rate, report_dict)
        return report_dict
