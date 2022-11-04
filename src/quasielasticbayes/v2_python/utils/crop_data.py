from numpy import ndarray
import numpy as np


def crop(x_data: ndarray, y_data: ndarray, e_data: ndarray,
         start_x: float, end_x: float) -> (ndarray, ndarray, ndarray):
    """
    Simple function for cropping x, y, e data
    """
    start_index = np.searchsorted(x_data, start_x)
    end_index = np.searchsorted(x_data, end_x)

    x_crop = x_data[start_index:end_index]
    y_crop = y_data[start_index:end_index]
    e_crop = None
    if e_data is not None:
        e_crop = e_data[start_index:end_index]
    return x_crop, y_crop, e_crop
