from quickBayes.functions.convolution import (
        ConvolutionWithResolution)
from quickBayes.functions.base import BaseFitFunction
from quickBayes.functions.delta import Delta
from numpy import ndarray
from typing import Dict, List
from abc import abstractmethod
import copy


"""
There are no direct tests for this class.
This is because all of the functionality
is tested by the classes that inherit
from it
"""


class QEFunction(BaseFitFunction):
    def __init__(self, bg_function: BaseFitFunction, elastic_peak: bool,
                 r_x: ndarray, r_y: ndarray, start_x: float, end_x: float):
        """
        Create a quasi elastic fitting function

        ASSUMPTIONS:
        - The added 'peak' is always the same (e.g. Lorentzian or stretch exp)
        - The 2nd param for the added function is the 'peak centre'
        - Will convolve with the resolution function
        - The fit function has more 'peaks' than the number in a report
        - The elastic peak can be represented by a 'delta' function

        :param bg_function: background fitting function
        :param elastic_peak: if to include an elastic peak (True/False)
        :param res_x: x values for resolution function
        :param res_y: y values for resolution function
        :param start_x: the start of the fitting range
        :param end_x: the end of the fitting range
        """
        self._N_peaks = 0
        self.BG = copy.copy(bg_function)
        self.BG.add_to_prefix(self.prefix + 'f1')
        self.conv = ConvolutionWithResolution(r_x, r_y, start_x,
                                              end_x, self.prefix + 'f2')
        self.delta = elastic_peak
        if elastic_peak:
            delta = Delta('')
            self.conv.add_function(delta)

        self._N_params = self.BG.N_params + self.conv.N_params,
        self._prefix = self.prefix

    def update_x_range(self, new_x: ndarray) -> None:
        """
        The sampling of the resolution function can make
        a big difference to the quality of the function.
        This is because of sampling issues. To solve
        this problem a user can update the x (and y)
        ranges using this method.
        :param new_x: the new x range (y is interpolated)
        """
        self.conv.update_x_range(new_x)

    @property
    def N_params(self) -> int:
        """
        :return the number of parameters in the function
        """
        # subtract 1 to share the peak position with delta
        peak_correction = self._N_peaks
        if not self.delta and self._N_peaks > 0:
            peak_correction -= 1
        return self.BG.N_params + self.conv.N_params - 1*peak_correction

    @property
    def N_peaks(self) -> int:
        """
        :return the number of extra functions (e.g. lorentzians, stretch exp)
        """
        return self._N_peaks

    @property
    def prefix(self) -> str:
        """
        :return the label for the number of peaks
        """
        return str(f'N{self.N_peaks}:')

    def _update_prefixes(self) -> None:
        """
        Method for updaing the prefixes for new number of peaks
        """
        self.BG.update_prefix(self.prefix)
        self.conv.update_prefix(self.prefix)

    def add_single_function(self, func: BaseFitFunction) -> None:
        """
        :param func: the function (e.g. lorentzian) to add
        Add a single Lorentzian function to the qldata function
        """
        self._N_peaks += 1
        self.conv.add_function(func)
        # update the labels/prefixes
        self._update_prefixes()

    def _add_params(self, offset: int, x0: float,
                    args: List[float]) -> List[float]:
        """
        :param offset: the (index) offset for the function
        :param x0: the peak centre
        :param args: the argument list
        :return the extended (with repeats) parameters to add
        """
        return []

    def _get_params(self, args: List[float]) -> List[float]:
        """
        For fitting we need to tie the peak centers for the
        delta and functions. This function creates the
        extended parameter list (with repeats).
        :param args: the arguments to the function (no repeats)
        :return the arguments with repeats for the peak centers
        in the correct places
        """
        params = []
        N_BG_params = self.BG.N_params
        x0 = 0
        N_f0 = 0
        offset = 0
        if len(self.conv._funcs) > 0:
            N_f0 = self.conv._funcs[0].N_params
            params = [*args[N_BG_params:N_BG_params + N_f0]]
            x0 = args[N_BG_params + 1]  # same position for both lor and delta
            if not self.delta:
                # if not elastic, already done first peak
                offset = 1

        for j in range(self._N_peaks - offset):
            params += self._add_params(N_BG_params + N_f0 + j*2, x0, args)
        return params

    def __call__(self, x: ndarray, *args: float) -> ndarray:
        """
        Implement the function evaluation.
        Need to follow the expected
        form for scipy
        :param x: x values for function evaluation
        :param args: args for functions
        :return y values for the function evaluation
        """
        N_BG_params = self.BG.N_params
        result = self.BG(x, *args[:N_BG_params])

        params = self._get_params(args)
        result += self.conv(x, *params)
        return result

    def _get_func_from_report(self, args: List[float]) -> List[float]:
        return args

    def read_from_report(self, report_dict: Dict[str, List[float]],
                         N: int, index: int = 0) -> List[float]:
        """
        Read the parameters from the results dict
        :param report_dict: the dict of results
        :param N: the number of peaks
        :param index: the index to get results from
        :return the parameters
        """
        if N > self._N_peaks:
            raise ValueError("Too many peaks selected")
        N_peaks = self._N_peaks
        # set for number of peaks
        self._N_peaks = N
        self._update_prefixes()
        # get parameters
        params = self.BG.read_from_report(report_dict, index)
        num_funcs = N
        if self.delta:
            num_funcs += 1

        if num_funcs > 0:
            # always want the first member in full
            params += self.conv._funcs[0].read_from_report(report_dict, index)

        for k in range(1, num_funcs):
            tmp = self.conv._funcs[k].read_from_report(report_dict,
                                                       index)
            params += self._get_func_from_report(tmp)

        # reset the number of peaks
        self._N_peaks = N_peaks
        self._update_prefixes()
        return params

    def report(self, report_dict: Dict[str, List[float]],
               *args: float) -> Dict[str, List[float]]:
        """
        Reports the results
        :param report_dict: dict of results
        :param args: args for functions
        :returns updated results dict
        """
        N = self.N_params
        if len(args) != N:
            raise ValueError(f"Expected {N} args, got {len(args)}")
        report_dict = self.BG.report(report_dict, *args[:self.BG.N_params])

        params = self._get_params(args)
        report_dict = self.conv.report(report_dict, *params)
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
        if len(errors) != self.N_params or len(params) != self.N_params:
            raise ValueError(f"Expected {self.N_params} args, "
                             f"got {len(params)} and {len(errors)}")
        report_dict = self.BG.report_errors(report_dict,
                                            errors[:self.BG.N_params],
                                            params[:self.BG.N_params])

        params = self._get_params(params)
        errors = self._get_params(errors)
        report_dict = self.conv.report_errors(report_dict, errors, params)
        return report_dict

    @abstractmethod
    def update_first_values(self, update: List[float],
                            guess: List[float], index=-1) -> None:
        """
        Update the values due to ties.
        This will depend on the function
        :param update: the values to update due to ties
        :param guess: the guess for the function
        :param index: the index of the function
        """
        raise NotImplementedError()

    def _func_guess(self, full_guess: List[float]) -> List[float]:
        """
        Extracts the guess for the function
        :param full_guess: the full list of guess parameters
        :return the reduced guess parameters (no repeats for peak centre)
        """
        return full_guess

    def get_guess(self) -> List[float]:
        """
        Get the intial guess values.
        This takes into account the tied
        parameters.
        :return a list of guess parameters for the fit
        """
        # need to copy to prevent incrementing guess by calling this function
        guess = copy.copy(self.BG.get_guess())
        # want to reduce the guess to remove tied parameters
        if len(self.conv._funcs) > 0:
            guess += copy.copy(self.conv._funcs[0].get_guess())
            for j in range(1, len(self.conv._funcs)):
                full_guess = copy.copy(self.conv._funcs[j].get_guess())
                guess += self._func_guess(full_guess)
        return guess

    def get_func_guess(self, index=-1):
        if self.N_peaks == 0:
            return
        delta_offset = 1 if self.delta and index != -1 else 0
        guess = self.conv._funcs[index + delta_offset].get_guess()
        tie_values = self.conv._funcs[0].get_guess()
        guess = self.update_first_values(guess, tie_values)
        return guess

    def set_guess(self, guess: List[float]) -> None:
        """
        Override the default guess setting.
        :param guess: the guess
        """
        raise RuntimeError("set_guess is not available for this "
                           "function. Please use: \n"
                           "- set_BG_guess \n"
                           "- set_delta_guess \n"
                           "- set_func_guess")

    def set_BG_guess(self, guess: List[float]) -> None:
        """
        Set the background guess values.
        :param guess: the guess for the BG
        """
        self.BG.set_guess(guess)

    def set_delta_guess(self, guess: List[float]) -> None:
        """
        Set the delta guess values.
        :param guess: the guess for the delta function
        """
        if self.delta:
            self.conv._funcs[0].set_guess(guess)

    def set_func_guess(self, guess: List[float], index=-1) -> None:
        """
        Set the  guess values.
        :param guess: the guess for the function
        :param index: the index of the function
        """
        # no func
        if self._N_peaks == 0:
            return
        delta_offset = 1 if self.delta and index != -1 else 0

        self.conv._funcs[index + delta_offset].set_guess(guess)
        to_update = self.conv._funcs[0].get_guess()
        to_update = self.update_first_values(to_update, guess)
        self.conv._funcs[0].set_guess(to_update)

    def get_bounds(self) -> (List[float], List[float]):
        """
        Gets the bounds for the parameters.
        :return a list of the lower and upper bounds
        """
        # need to copy to prevent incrementing bounds by calling this function
        bounds = self.BG.get_bounds()
        lower = copy.copy(bounds[0])
        upper = copy.copy(bounds[1])

        if len(self.conv._funcs) > 0:
            # want to reduce the guess to remove tied parameters
            func = self.conv._funcs[0]
            bounds = func.get_bounds()
            lower += copy.copy(bounds[0])
            upper += copy.copy(bounds[1])

            for j in range(1, len(self.conv._funcs)):
                f = self.conv._funcs[j]
                bounds = f.get_bounds()
                tmp = bounds[0]
                lower += copy.copy(self._func_guess(tmp))
                tmp = bounds[1]
                upper += copy.copy(self._func_guess(tmp))

        return lower, upper

    def set_bounds(self, lower: List[float], upper: List[float]) -> None:
        """
        Override the default bounds setting.
        :param lower: the lower bound values
        :param upper: the upper bound values
        """
        raise RuntimeError("set_guess is not available for this "
                           "function. Please use: \n"
                           "- set_BG_bounds \n"
                           "- set_delta_bounds \n"
                           "- set_func_bounds")

    def set_BG_bounds(self, lower: List[float], upper: List[float]) -> None:
        """
        Set the background bounds values.
        :param lower: the lower bound values
        :param upper: the upper bound values
        """
        self.BG.set_bounds(lower, upper)

    def set_delta_bounds(self, lower: List[float], upper: List[float]) -> None:
        """
        Set the delta bounds values.
        :param lower: the lower bound values
        :param upper: the upper bound values
        """
        if self.delta:
            self.conv._funcs[0].set_bounds(lower, upper)

    def set_func_bounds(self, lower: List[float],
                        upper: List[float], index=-1) -> None:
        """
        Set the bounds values.
        :param lower: the lower bound values
        :param upper: the upper bound values
        :param index: the index of the function
        """
        if self._N_peaks == 0:
            return
        delta_offset = 1 if self.delta and index != -1 else 0

        self.conv._funcs[index + delta_offset].set_bounds(lower, upper)
        update_lower, update_upper = self.conv._funcs[0].get_bounds()
        update_lower = self.update_first_values(update_lower, lower)
        update_upper = self.update_first_values(update_upper, upper)
        self.conv._funcs[0].set_bounds(update_lower, update_upper)
