from quasielasticbayes.v2.base import BaseFitFunction
from numpy import ndarray, pi
from typing import Dict, List


class Lorentzian(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Create a Lorentzian function
        :param prefix: the prefix for the parameters
        """
        super().__init__(3, prefix)

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
        report_dict = self._add_to_report(f"{self._prefix}Amplitude",
                                          a, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Peak Centre",
                                          x0, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Gamma",
                                          Gamma, report_dict)
        return report_dict
