from quickBayes.functions.base import BaseFitFunction
from quickBayes.functions.composite import CompositeFunction
from quickBayes.utils.crop_data import crop
from quickBayes.utils.spline import spline
from numpy import ndarray
from scipy import signal
import copy


class ConvolutionWithResolution(CompositeFunction):
    def __init__(self, res_x: ndarray, res_y: ndarray,
                 start_x: float, end_x: float, prefix: str = ''):
        """
        Creates a convolution with a tabulated resolution function.
        Can add functions to be convoluted with the resolution.
        :param res_x: x values for resolution function
        :param res_y: y values for resolution function
        :param start_x: the start of the fitting range
        :param end_x: the end of the fitting range
        :param prefix: prefix for fitting function
        """
        super().__init__(prefix)
        self._rx = copy.deepcopy(res_x)
        self._ry = copy.deepcopy(res_y)
        self._rx, self._ry, _ = crop(self._rx, self._ry, None, start_x, end_x)
        # this is to normalise the kernal to get correct amplitudes
        self._ry /= sum(self._ry)

    def update_x_range(self, new_x: ndarray) -> None:
        """
        The sampling of the resolution function can make
        a big difference to the quality of the function.
        This is because of sampling issues. To solve
        this problem a user can update the x (and y)
        ranges using this method.
        :param new_x: the new x range (y is interpolated)
        """
        self._ry = spline(self._rx, self._ry, new_x)
        self._ry /= sum(self._ry)
        self._rx = new_x

    def update_prefix(self, new: str) -> None:
        """
        Update the begining of the prefixes
        :param new: the new part of the prefix
        """
        for j in range(len(self._funcs)):
            self._funcs[j].update_prefix(new)

    def add_function(self, func: BaseFitFunction) -> None:
        """
        Adds a function to the convolution
        :param func: the function to add
        """
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
        :param x: x range to calculate function over
        :param args: the arguments for the convolution function
        :return y values for the convolution
        """
        result = super().__call__(x, *args)
        # assume rx and x are the same
        return signal.convolve(result, self._ry, mode='same')
