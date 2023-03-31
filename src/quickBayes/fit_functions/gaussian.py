from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray, pi
import numpy as np
from typing import Dict, List


class Gaussian(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Create a gaussian function
        :param prefix: prefix for parameter reporting
        """
        super().__init__(3, prefix, [1., 0, 0.1],
                         [0., -1., 0.], [np.inf, 1., np.inf])

    @property
    def amplitude(self) -> str:
        """
        :retrun string for amplitude name
        """
        return str(f"{self._prefix}Amplitude")

    @property
    def mean(self) -> str:
        """
        :return string for name of mean
        """
        return str(f"{self._prefix}Mean")

    @property
    def sigma(self) -> str:
        """
        :return string for the name of sigma
        """
        return str(f"{self._prefix}Sigma")

    def __call__(self, x: ndarray, amplitude: float, x0: float,
                 sigma: float) -> ndarray:
        """
        Implement a gaussian.
        Need to follow the expected
        form for scipy
        :param x: x values for the function evaluation
        :param amplitude: amplitude of gaussian
        :param x0: the mean value of the gaussian
        :param sigma: the sigma value of the gaussian
        :return y values for the gaussian
        """
        pre_factor = amplitude/(sigma*np.sqrt(2.*pi))
        return pre_factor*np.exp(-pow(x-x0, 2)/(2.*sigma*sigma))

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.mean, index),
                self._read_report(report_dict, self.sigma, index)]

    def report(self, report_dict: Dict[str, List[float]], a: float,
               x0: float, sigma: float) -> Dict[str, List[float]]:
        """
        Report the results
        :param report_dict: dict of results
        :param a: amplitude of gaussian
        :param x0: the mean value of the gaussian
        :param sigma: the sigma value of the gaussian
        returns the updated results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          a, report_dict)
        report_dict = self._add_to_report(self.mean,
                                          x0, report_dict)
        report_dict = self._add_to_report(self.sigma,
                                          sigma, report_dict)
        return report_dict
