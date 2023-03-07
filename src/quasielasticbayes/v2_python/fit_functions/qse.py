from quasielasticbayes.v2.functions.base import BaseFitFunction
from quasielasticbayes.v2.functions.qe_function import QEFunction
from quasielasticbayes.v2.functions.SE import StretchExp
from numpy import ndarray
import copy
from typing import List


class QSEFunction(QEFunction):
    def __init__(self, bg_function: BaseFitFunction, elastic_peak: bool,
                 r_x: ndarray, r_y: ndarray, start_x: float, end_x: float):
        """
        Create a quasi elastic fitting function using stretched exp
        :param bg_function: background fitting function
        :param elastic_peak: if to include an elastic peak (True/False)
        :param res_x: x values for resolution function
        :param res_y: y values for resolution function
        :param start_x: the start of the fitting range
        :param end_x: the end of the fitting range
        """
        super().__init__(bg_function, elastic_peak, r_x,
                         r_y, start_x, end_x)

    def add_single_SE(self) -> None:
        """
        Add a single Lorentzian function to the qldata function
        """
        se = StretchExp()
        self.add_single_function(se)

    def _add_params(self, offset: int, x0: float,
                    args: List[float]) -> List[float]:

        # adds amplitude, peak centre and FWHM
        return [args[offset], x0, args[offset + 1], args[offset + 2]]

    def _get_func_from_report(self, args: List[float]) -> List[float]:
        # skip the peak centre
        return [args[0], args[2], args[3]]

    def _func_guess(self, full_guess: List[float]) -> List[float]:
        """
        Get the intial guess values.
        This takes into account the tied
        parameters.
        :return a list of guess parameters for the fit
        """
        # skip peak centre
        return [full_guess[0], full_guess[2], full_guess[3]]

    def get_guess(self) -> List[float]:
        """
        Gets the guess for the fit params
        :result a list of initial values for fit
        """
        guess = copy.copy(self.BG.get_guess())

        # want to reduce the guess to remove tied paramaters
        if len(self.conv._funcs) > 0 and self.delta:
            guess += copy.copy(self.conv._funcs[0].get_guess())
        elif len(self.conv._funcs) > 0:
            guess += copy.copy(self.conv._funcs[0].get_guess())
        if len(self.conv._funcs) > 1:
            for j in range(1, len(self.conv._funcs)):
                full_guess = copy.copy(self.conv._funcs[j].get_guess())
                guess += self._func_guess(full_guess)
        return guess

    def update_first_values(self, to_update: List[float],
                            guess: List[float]) -> List[float]:
        """
        Method for copying the updated values into the first function
        in the convolution (this determines the value in evaluation.
        :param to_update: the values to update (due to ties)
        :param guess: the new guess values for the function being changed
        :return the updated list
        """
        to_update[1] = guess[1]
        return to_update
