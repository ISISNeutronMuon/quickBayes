from numpy import ndarray
from scipy.interpolate import interp1d


def spline(x_data: ndarray, y_data: ndarray,
           new_x_values: ndarray) -> ndarray:
    """
    Simple function for interpolating and extrapolating
    (to zero) data.
    :param x_data: the origianl x data
    :param y_data: the original y data
    :param new_x_values: the new x data
    :return the new y values
    """
    func = interp1d(x_data, y_data, bounds_error=False,
                    fill_value=0., kind='cubic')
    return func(new_x_values)
