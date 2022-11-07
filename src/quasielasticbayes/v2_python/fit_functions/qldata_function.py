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
        self.BG = bg_function
        self.BG.add_to_prefix('f1')
        self.conv = ConvolutionWithResolution(r_x, r_y, start_x, end_x, 'f2')
        self._N_peaks = 0
        self.elastic = elastic_peak

        self.delta = elastic_peak
        if elastic_peak:
            delta = Delta('')
            self.conv.add_function(delta)
        super().__init__(0, '')

    @property
    def N_params(self):
        """
        :return the number of parameters in the function
        """
        # subtract 1 to share the peak position with delta
        return self.BG.N_params + self.conv.N_params - 1

    def add_single_lorentzian(self):
        """
        Add a single Lorentzian function to the qldata function
        """
        lor = Lorentzian(self._prefix)
        self.conv.add_function(lor)
        self._N_peaks += 1

    def _get_params(self, args: List[float]) -> List[float]:
        # tie the peak centers
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
                qe_amp = args[BG_N_params + N_e + j*N_qe]
                EISF = e_amp/(e_amp + qe_amp)
                report_dict = self._add_to_report(self.conv._funcs[k]._prefix
                                                  + "EISF",
                                                  EISF, report_dict)
        return report_dict
