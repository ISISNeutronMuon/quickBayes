from numpy import ndarray
import numpy as np


def crop(x_data: ndarray, y_data: ndarray, e_data: ndarray,
         start_x: float, end_x: float) -> (ndarray, ndarray, ndarray):
    """
    Simple function for cropping x, y, e data
    :param x_data: x data to crop
    :param y_data: y data to crop
    :param e_data: error data to crop, does nothing if None
    :param start_x: x value to start cropping from
    :param end_x: x value to stop cropping at
    :return cropped values for x, y, e (None if input is None)
    """
    start_index = np.searchsorted(x_data, start_x)
    end_index = np.searchsorted(x_data, end_x)

    x_crop = x_data[start_index:end_index]
    y_crop = y_data[start_index:end_index]
    e_crop = None
    if e_data is not None:
        e_crop = e_data[start_index:end_index]
    return x_crop, y_crop, e_crop
