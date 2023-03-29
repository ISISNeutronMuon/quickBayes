from quickBayes.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List
from scipy.fftpack import fft, fftfreq
from scipy.special import gamma
from scipy import constants


"""
This code is taken from the open source code
Mantid.
https://github.com/mantidproject/mantid/blob/main/Framework/PythonInterface/plugins/functions/StretchedExpFTHelper.py
"""


PLANCK_CONSTANT = constants.Planck / constants.e * 1.e15  # meV*psec


def function1Dcommon(xvals: ndarray, tau: float, beta: float,
                     refine_factor=16,) -> (ndarray, ndarray):
    """
    Fourier transform of the symmetrized stretched exponential
    :param function: instance of StretchedExpFT
    :param xvals: energy domain
    :param tau: relaxation time
    :param beta: stretching exponenet
    :param refine_factor: divide the natural energy width by this value
    :return: energies, and function values
    """
    N = len(xvals)

    # energy spacing. Assumed xvals is a single-segment grid
    # of increasing energy values
    dE = (xvals[-1] - xvals[0]) / (refine_factor * (N - 1))
    E_range = 2 * max(abs(xvals))

    dt = 0.5 * PLANCK_CONSTANT / E_range  # spacing in time
    tmax = PLANCK_CONSTANT / dE  # maximum reciprocal time
    # round to an upper power of two
    nt = 2 ** (1 + int(np.log(tmax / dt) / np.log(2)))
    sampled_times = dt * np.arange(-nt, nt)

    decay = np.exp(-(np.abs(sampled_times) / tau)**beta)

    """
    The Fourier transform introduces an extra factor exp(i*pi*E/dE),
    which amounts to alternating sign every time E increases by dE,
    the energy bin width. Thus, we take the absolute value
    """
    fourier = np.abs(fft(decay).real)  # notice the reverse of decay array

    fourier /= fourier[0]  # set maximum to unity
    # Normalize the integral in energies to unity
    fourier *= 2*tau*gamma(1./beta) / (beta*PLANCK_CONSTANT)

    # symmetrize to negative energies
    fourier = np.concatenate(
        [fourier[nt:], fourier[:nt]])  # increasing ordering

    # Find energy values corresponding to the fourier values
    energies = PLANCK_CONSTANT * fftfreq(2 * nt, d=dt)  # standard ordering
    energies = np.concatenate(
        [energies[nt:], energies[:nt]])  # increasing ordering
    return energies, fourier


class StretchExp(BaseFitFunction):
    def __init__(self, prefix: str = ''):
        """
        Create a stretched exponenetial function
        :param prefix: the prefix for the parameters
        """
        super().__init__(4, prefix,
                         [.1, 0.0, 6.582, 0.7],  # 6.582 -> FWHM = 0.2
                         [0., -1., 0, 0], [1., 1., 100., 1.])

    @property
    def amplitude(self) -> str:
        """
        :return the string for the amplitude
        """
        return str(f"{self._prefix}Amplitude")

    @property
    def x0(self) -> str:
        """
        :return the string for the peak centre
        """
        return str(f"{self._prefix}Peak Centre")

    @property
    def beta(self) -> str:
        """
        :return the string for the beta value
        """
        return str(f"{self._prefix}beta")

    @property
    def tau_str(self) -> str:
        """
        :return the string value for tau
        """
        return str(f"{self._prefix}tau")

    def __call__(self, x: ndarray,
                 amplitude: float, x0: float,
                 tau: float, beta: float) -> ndarray:
        """
        Implement the stretched exponential.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param amplitude: amplitude
        :param x0: the peak centre
        :param tau: relaxation time
        :param beta: stretching exponent
        :return y values for function evaluation
        """

        energies, fourier = function1Dcommon(x, tau, beta)
        return amplitude*np.interp(x - x0,
                                   energies, fourier)

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param index: the index to get results from
        :return the parameters
        """
        return [self._read_report(report_dict, self.amplitude, index),
                self._read_report(report_dict, self.x0, index),
                self._read_report(report_dict, self.tau_str, index),
                self._read_report(report_dict, self.beta, index)]

    @staticmethod
    def FWHM(tau: float) -> float:
        """
        Method to convert tau to FWHM
        :param tau: tau parameter
        :return FWHM
        """
        return 2.*PLANCK_CONSTANT/(2.*np.pi*tau)

    @staticmethod
    def tau(FWHM: float) -> float:
        """
        Method to get tau from FWHM
        Used for estimation
        :param FWHM: full width half maximum
        :return tau parameter
        """
        return 2.*PLANCK_CONSTANT/(2.*np.pi*FWHM)

    def report(self, report_dict: Dict[str, List[float]], a: float, x0: float,
               tau: float, beta: float) -> Dict[str, List[float]]:
        """
        report parameters
        :param report_dic: dict of results
        :param a: amplitude
        :param x0: the peak centre
        :param tau: relaxation time
        :param beta: stretching exponent
        :return update results dict
        """
        report_dict = self._add_to_report(self.amplitude,
                                          a, report_dict)
        report_dict = self._add_to_report(self.x0,
                                          x0, report_dict)
        report_dict = self._add_to_report(self.beta,
                                          beta, report_dict)
        report_dict = self._add_to_report(self.tau_str,
                                          tau, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}FWHM",
                                          self.FWHM(tau), report_dict)
        return report_dict

    def report_errors(self, report_dict: Dict[str, List[float]],
                      errors: ndarray,
                      params: ndarray) -> Dict[str, List[float]]:
        """
        report parameters
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
                                          errors[3], report_dict)
        report_dict = self._add_to_report(self.tau_str,
                                          errors[2], report_dict)
        """
        FWHM = 2 hbar /tau
        sigma_{FWHM} = 2 hbar sigma_tau / tau^2 = FWHM sigma_tau/tau
        """
        tmp = errors[2]/params[2]
        report_dict = self._add_to_report(f"{self._prefix}FWHM",
                                          self.FWHM(params[2])*tmp,
                                          report_dict)
        return report_dict

    def set_guess_FWHM(self, values: List[float]) -> List[float]:
        """
        set the starting guess for a fit function
        :param values: the guess values. The 3rd value is the FWHM
        """
        self._check_length(values, 'guess')
        self._guess = [values[0], values[1],
                       self.tau(values[2]), values[3]]

    def set_guess(self, values: List[float]) -> List[float]:
        """
        set the starting guess for a fit function
        :param values: the guess values. The 3rd value is the FWHM
        """
        self._check_length(values, 'guess')
        self._guess = values
