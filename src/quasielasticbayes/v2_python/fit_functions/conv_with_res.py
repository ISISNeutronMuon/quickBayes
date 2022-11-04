from quasielasticbayes.v2.composite import CompositeFunction
from quasielasticbayes.v2.crop_data import crop
from numpy import ndarray
from scipy import signal


class ConvolutionWithResolution(CompositeFunction):
    def __init__(self, res_x: ndarray, res_y: ndarray,
                 start_x: float, end_x: float, prefix: str = ''):
        super().__init__(prefix)
        self._rx, self._ry, _ = crop(res_x, res_y, None, start_x, end_x)
        # this is to normalise the kernal to get correct amplitudes
        self._ry /= sum(self._ry)

    def add_function(self, func):
        # add prefix for convolution
        if self._prefix == '':
            self._prefix = 'f1'
        func.add_to_prefix(self._prefix)
        super().add_function(func)

    def __call__(self, x: ndarray, *args: float) -> ndarray:
        """
        Implement a convolution with a resolution function.
        Need to follow the expected
        form for scipy
        """
        result = super().__call__(x, *args)
        # assume rx and x are the same
        return signal.convolve(result, self._ry, mode='same')
