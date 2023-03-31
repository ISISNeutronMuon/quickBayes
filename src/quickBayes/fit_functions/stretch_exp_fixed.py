from quickBayes.functions.SE import StretchExp
from numpy import ndarray
from typing import Dict, List


class StretchExpWithFixes(StretchExp):
    def __init__(self, FWHM: float = 0.2, beta: float = 0.8, prefix: str = ''):
        """
        Create a stretched exponential function with 2 fixed parameters.
        :param FWHM: full width half max value for the fix
        :param beta: the beta value for the fix
        :param prefix: the prefix for the parameters
        """
        self._func = StretchExp()
        self.set_beta(beta)
        self.set_FWHM(FWHM)
        super().__init__(prefix)
        # change stuff for 2 free parameters
        self._N_params = 2
        self._guess = self._guess[0:2]
        self._lower = self._lower[0:2]
        self._upper = self._upper[0:2]

    def set_FWHM(self, FWHM: float) -> None:
        """
        Update the FWHM fix value
        :param FWHM: full width half max for fix
        """
        self._tau = self.tau(FWHM)

    def set_beta(self, beta: float) -> None:
        """
        Update the beta fix value
        :param beta: the beta value for fix
        """
        self._beta = beta

    @property
    def get_tau(self) -> float:
        """
        Get the tau value being used.
        tau is related to the FWHM.
        :return the tau value used in fix
        """
        return self._tau

    @property
    def get_beta(self) -> float:
        """
        Gets the beta value being used
        :return beta value for fix
        """
        return self._beta

    def __call__(self, x: ndarray,
                 amplitude: float, x0: float) -> ndarray:
        """
        Implement the stretched exponential.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param amplitude: amplitude
        :param x0: the peak centre
        :return y values for function evaluation
        """
        return super().__call__(x, amplitude, x0, self.get_tau, self.get_beta)

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        and sets beta and tau
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        self._tau = self._read_report(report_dict, self.tau_str, index)
        self.set_beta(self._read_report(report_dict, self.beta, index))

        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.x0, index)]

    def set_guess_FWHM(self, value: List[float]) -> None:
        """
        This is an inherited function that will not do anything
        :param value: value to set
        """
        raise RuntimeError("this method does not work with fix")

    def report(self, report_dict: Dict[str, List[float]],
               a: float, x0: float) -> Dict[str, List[float]]:
        """
        report parameters, including the fixed ones.
        :param report_dic: dict of results
        :param a: amplitude
        :param x0: the peak centre
        :return update results dict
        """
        return super().report(report_dict, a, x0,
                              self.get_tau, self.get_beta)

    def report_errors(self, report_dict: Dict[str, List[float]],
                      errors: ndarray,
                      params: ndarray) -> Dict[str, List[float]]:
        """
        report parameters. The errors are 0 for the
        fixed parameters.
        :param report_dic: dict of parameter errors
        :param errors: the errors for the fit parameters
        :param params: the fit parameters
        :return update results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          errors[0], report_dict)
        report_dict = self._add_to_report(self.x0,
                                          errors[1], report_dict)
        report_dict = self._add_to_report(self.beta,
                                          0.0, report_dict)
        report_dict = self._add_to_report(self.tau_str,
                                          0.0, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}FWHM",
                                          0,
                                          report_dict)
        return report_dict
