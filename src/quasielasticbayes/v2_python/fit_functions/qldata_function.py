from quasielasticbayes.v2.convolution import ConvolutionWithResolution
from quasielasticbayes.v2.base import BaseFitFunction
from quasielasticbayes.v2.delta import Delta
from quasielasticbayes.v2.lorentz import Lorentzian
from numpy import ndarray
from typing import Dict, List


class QlDataFunction(BaseFitFunction):
    def __init__(self, bg_function: BaseFitFunction, elastic_peak: bool,
                 r_x: ndarray, r_y: ndarray, start_x: float, end_x: float):
        self.BG = bg_function
        self.BG.add_to_prefix('f1')
        self.conv = ConvolutionWithResolution(r_x, r_y, start_x, end_x, 'f2')
        if elastic_peak:
            delta = Delta('')
            self.conv.add_function(delta)
        super().__init__(0, '')

    @property
    def N_params(self):
        return self.BG.N_params + self.conv.N_params

    def add_single_lorentzian(self):
        lor = Lorentzian(self._prefix)
        self.conv.add_function(lor)

    def __call__(self, x: ndarray, *args: float) -> ndarray:
        """
        Implement the function evaluation.
        Need to follow the expected
        form for scipy
        """
        result = self.BG(x, *args[:self.BG.N_params])
        result += self.conv(x, *args[self.BG.N_params:])
        return result

    def report(self, report_dict: Dict[str, List[float]],
               *args: float) -> Dict[str, List[float]]:
        """
        returns the fit parameters as a dict
        """
        if len(args) != self.N_params:
            raise ValueError(f"Expected {self.N_params} args, got {len(args)}")
        report_dict = self.BG.report(report_dict, *args[:self.BG.N_params])
        report_dict = self.conv.report(report_dict, *args[self.BG.N_params:])
        return report_dict
