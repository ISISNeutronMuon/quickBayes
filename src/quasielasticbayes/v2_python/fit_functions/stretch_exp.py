from quasielasticbayes.v2.functions.base import BaseFitFunction
from numpy import ndarray
import numpy as np
from typing import Dict, List
from scipy.fftpack import fft, fftfreq
from scipy.special import gamma
from scipy import constants


"""
This code is taken from the open source code
Mantid.
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
        super().__init__(4, prefix)

    def __call__(self, x: ndarray, x0: float,
                 amplitude: float,
                 tau: float, beta: float) -> ndarray:
        """
        Implement the stretched exponential.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param x0: the peak centre
        :param amplitude: amplitude
        :param tau: relaxation time
        :param beta: stretching exponent
        :return y values for function evaluation
        """

        energies, fourier = function1Dcommon(x, tau, beta)
        return amplitude*np.interp(x - x0,
                                   energies, fourier)

    def report(self, report_dict: Dict[str, List[float]], x0: float, a: float,
               tau, beta) -> Dict[str, List[float]]:
        """
        report parameters
        :param report_dic: dict of results
        :param x0: the peak centre
        :param a: amplitude
        :param tau: relaxation time
        :param beta: stretching exponent
        :return update results dict
        """
        report_dict = self._add_to_report(f"{self._prefix}Amplitude",
                                          a, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}Peak Centre",
                                          x0, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}beta",
                                          beta, report_dict)
        report_dict = self._add_to_report(f"{self._prefix}tau",
                                          tau, report_dict)
        return report_dict

    def get_guess(self) -> List[float]:
        """
        Get the starting guess for a fit function
        :return the initial guess
        """
        return [.0, 0.01, 0.01, 0.01]

    def get_bounds(self) -> (List[float], List[float]):
        """
        Get the fitting bounds
        :return lists for lower and upper bounds
        """
        return [-1, 0., 0, 0], [1., 1., 100., 1.]
