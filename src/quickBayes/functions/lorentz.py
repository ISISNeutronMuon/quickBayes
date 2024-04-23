from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray, pi
from typing import Dict, List


class Lorentzian(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Create a Lorentzian function
        :param prefix: the prefix for the parameters
        """
        super().__init__(3, prefix, [0.01, 0.0, 0.02],
                         [0., -1., 1.e-6], [1., 1., 1.])

    @property
    def amplitude(self) -> str:
        """
        :return string for the amplitude
        """
        return str(f"{self._prefix}Amplitude")

    @property
    def centre(self) -> str:
        """
        :return the string for the peak centre
        """
        return str(f"{self._prefix}Peak Centre")

    @property
    def Gamma(self) -> str:
        """
        :return the string for Gamma
        """
        return str(f"{self._prefix}Gamma")

    def __call__(self, x: ndarray, amplitude: float,
                 x0: float, Gamma: float) -> ndarray:
        """
        Implement the Lorentzian.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param amplitude: amplitude of the lorentzian
        :param x0: the peak centre
        :param Gamma: half width at half maxima (HWHM)
        :return y values for function evaluation
        """
        G = Gamma/2.
        return amplitude*G/(pi*(pow(x-x0, 2)+pow(G, 2)))

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.centre, index),
                self._read_report(report_dict, self.Gamma, index)]

    def report(self, report_dict: Dict[str, List[float]], a: float,
               x0: float, Gamma: float) -> Dict[str, List[float]]:
        """
        report parameters
        :param report_dic: dict of results
        :param a: amplitude of the lorentzian
        :param x0: the peak centre
        :param Gamma: half width at half maxima (HWHM)
        :return update results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          a, report_dict)
        report_dict = self._add_to_report(self.centre,
                                          x0, report_dict)
        report_dict = self._add_to_report(self.Gamma,
                                          Gamma, report_dict)
        return report_dict
