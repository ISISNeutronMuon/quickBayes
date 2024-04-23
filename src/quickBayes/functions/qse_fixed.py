from quickBayes.functions.SE_fix import StretchExpWithFixes
from quickBayes.functions.base import BaseFitFunction
from quickBayes.functions.qse_function import QSEFunction
from numpy import ndarray
from typing import List


class QSEFixFunction(QSEFunction):

    def __init__(self, bg_function: BaseFitFunction, elastic_peak: bool,
                 r_x: ndarray, r_y: ndarray, start_x: float, end_x: float):
        """
        Create a quasi elastic fitting function using fixed stretched exp
        :param bg_function: background fitting function
        :param elastic_peak: if to include an elastic peak (True/False)
        :param res_x: x values for resolution function
        :param res_y: y values for resolution function
        :param start_x: the start of the fitting range
        :param end_x: the end of the fitting range
        """
        self._se = []
        super().__init__(bg_function, elastic_peak, r_x,
                         r_y, start_x, end_x)

    def add_single_SE(self) -> None:
        """
        adds a single stretched exp with fixes
        """
        self._se.append(StretchExpWithFixes())
        self.add_single_function(self._se[-1])

    @staticmethod
    def _get_func_from_report(args: List[float]) -> List[float]:
        """
        extracts the relevant info from report (excluding fixes)
        :param args: the full list of arguments.
        :return the arguments, skipping the peak centre
        """
        return [args[0]]

    @staticmethod
    def _func_guess(full_guess: List[float]) -> List[float]:
        """
        Get the initial guess values.
        This takes into account the tied
        parameters.
        :return a list of guess parameters for the fit
        """
        # skip peak centre
        return [full_guess[0]]

    def set_beta(self, beta: float, index: int = -1) -> None:
        """
        sets the beta value for the fix
        :param beta: the beta value to fix to
        :param index: the index of the se
        """
        if not self._se:
            return
        self._se[index].set_beta(beta)

    def set_FWHM(self, FWHM, index: int = -1) -> None:
        """
        sets the FWHM value for the fix
        :param FWHM: the FWHM value to fix to
        :param index: the index of the se
        """
        if not self._se:
            return
        self._se[index].set_FWHM(FWHM)

    @staticmethod
    def _add_params(offset: int, x0: float,
                    args: List[float]) -> List[float]:
        """
        Gets the amplitude and peak centre
        :param offset: the offset in the index
        :param x0: the peak centre
        :param args: the arguments
        :return a list of parameters for the SE
        """
        # adds amplitude, peak centre
        return [args[offset], x0]

    def set_func_guess_FWHM(self, guess: List[float], index=-1) -> None:
        """
        Set the  guess values.
        :param guess: the guess for the function, assume
        FWHM is in index 2
        :param index: the index of the function
        """
        if self.N_peaks == 0:
            return

        self.set_FWHM(guess[2], index)

        super().set_func_guess(guess[0:2], index)
