from quasielasticbayes.v2.functions.convolution import (
        ConvolutionWithResolution)
from quasielasticbayes.v2.functions.base import BaseFitFunction
from quasielasticbayes.v2.functions.delta import Delta
from quasielasticbayes.v2.functions.lorentz import Lorentzian
from numpy import ndarray
from typing import Dict, List


class QlDataFunction(BaseFitFunction):
    def __init__(self, bg_function: BaseFitFunction, elastic_peak: bool,
                 r_x: ndarray, r_y: ndarray, start_x: float, end_x: float):
        """
        Create a quasi elastic fitting function using lorentzians
        :param bg_function: background fitting function
        :param elastic_peak: if to include an elastic peak (True/False)
        :param res_x: x values for resolution function
        :param res_y: y values for resolution function
        :param start_x: the start of the fitting range
        :param end_x: the end of the fitting range
        """
        self._N_peaks = 0
        self.BG = bg_function
        self.BG.add_to_prefix(self.prefix + 'f1')
        self.conv = ConvolutionWithResolution(r_x, r_y, start_x,
                                              end_x, self.prefix + 'f2')
        self.elastic = elastic_peak

        self.delta = elastic_peak
        if elastic_peak:
            delta = Delta('')
            self.conv.add_function(delta)
        super().__init__(0, self.prefix)

    @property
    def N_params(self) -> int:
        """
        :return the number of parameters in the function
        """
        # subtract 1 to share the peak position with delta
        return self.BG.N_params + self.conv.N_params - 1*self._N_peaks

    @property
    def N_peaks(self) -> int:
        return self._N_peaks

    @property
    def prefix(self) -> str:
        return str(f'N{self.N_peaks}:')

    def add_single_lorentzian(self) -> None:
        """
        Add a single Lorentzian function to the qldata function
        """
        self._N_peaks += 1
        lor = Lorentzian()
        self.conv.add_function(lor)
        # update the labels/prefixes
        self.BG.update_prefix(self.prefix)
        self.conv.update_prefix(self.prefix)

    def _get_params(self, args: List[float]) -> List[float]:
        """
        For fitting we need to tie the peak centers for the
        delta and lorentzians. This function creates the
        extended parameter list (with repeats).
        :param args: the arguments to the function (no repeats)
        :return the arguments with repeats for the peak centers
        in the correct places
        """
        params = []
        k = 0
        N_BG_params = self.BG.N_params
        if self.elastic:
            params = [*args[N_BG_params:N_BG_params + 2]]  # 2 for delta
            k = 2
            x0 = args[N_BG_params + 1]
        elif self._N_peaks > 0:
            x0 = args[N_BG_params + 1]
        for j in range(self._N_peaks):
            params += [args[N_BG_params + k + j*2]]  # adds amp
            params += [x0]
            params += [args[N_BG_params + k + j*2 + 1]]  # adds FWHM
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
        # manually add EISF
        if self.delta and self.conv.N_params > 2:
            BG_N_params = self.BG.N_params
            e_amp = args[BG_N_params]
            N_e = self.conv._funcs[0].N_params
            for j in range(len(self.conv._funcs)-1):
                k = j + 1
                N_qe = self.conv._funcs[k].N_params
                # params has values for delta + lorentz
                qe_amp = params[j*N_qe + N_e]
                EISF = e_amp/(e_amp + qe_amp)
                report_dict = self._add_to_report(self.conv._funcs[k]._prefix
                                                  + "EISF",
                                                  EISF, report_dict)
        return report_dict

    def get_guess(self) -> List[float]:
        """
        Get the intial guess values.
        This takes into account the tied
        parameters.
        :return a list of guess parameters for the fit
        """
        guess = self.BG.get_guess()

        # want to reduce the guess to remove tied paramaters
        if len(self.conv._funcs) > 0:
            guess += self.conv._funcs[0].get_guess()
            for j in range(1, len(self.conv._funcs)):
                full_guess = self.conv._funcs[j].get_guess()
                guess += [full_guess[0], full_guess[2]]
        return guess

    def get_bounds(self) -> (List[float], List[float]):
        """
        Gets the bounds for the parameters.
        :return a list of the lower and upper bounds
        """
        bounds = self.BG.get_bounds()
        lower = bounds[0]
        upper = bounds[1]

        if len(self.conv._funcs) > 0:
            # want to reduce the guess to remove tied paramaters
            f = self.conv._funcs[0]
            bounds = f.get_bounds()
            lower += bounds[0]
            upper += bounds[1]

            for j in range(1, len(self.conv._funcs)):
                f = self.conv._funcs[j]
                bounds = f.get_bounds()
                tmp = bounds[0]
                lower += [tmp[0], tmp[2]]
                tmp = bounds[1]
                upper += [tmp[0], tmp[2]]

        return lower, upper
