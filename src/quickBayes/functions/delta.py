from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List


class Delta(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Strictly this is not a true delta function.
        Instead it is a top hat function, which
        in the limit of binwidth-> 0 is a delta
        :param prefix: prefix for the parameters
        """
        super().__init__(2, prefix, [1., 0.], [0., -1], [np.inf, 1.])

    @property
    def amplitude(self) -> str:
        """
        :return the string for the amplitude
        """
        return str(f"{self._prefix}Amplitude")

    @property
    def centre(self) -> str:
        """
        :return the string for peak centre
        """
        return str(f"{self._prefix}Centre")

    def __call__(self, x: ndarray, amplitude: float, x0: float) -> ndarray:
        """
        Implement the delta/top hat.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param amplitude: height of the top hat function
        :param x0: the position of the top hat
        :return y values for the function evaluation
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

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.centre, index)]

    def report(self, report_dict: Dict[str, List[float]],
               amplitude: float, x0: float) -> Dict[str, List[float]]:
        """
        report the results
        :param report_dict: dict of results
        :param amplitude: height of the top hat function
        :param x0: the position of the top hat
        :return updated results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          amplitude, report_dict)
        report_dict = self._add_to_report(self.centre,
                                          x0, report_dict)
        return report_dict
