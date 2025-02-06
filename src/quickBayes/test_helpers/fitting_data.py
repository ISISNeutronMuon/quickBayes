from numpy import ndarray
import numpy as np


"""
This provide the fitting data and functions
for the fitting test template.
"""


def func(x_data: ndarray) -> None:
    """
    :param x_data: the x data
    """
    y = 0.1 + 0.9*x_data
    y[-1] = 0.99*y[-1]
    return y


def basic_data() -> (ndarray, ndarray, ndarray):
    """
    Gets simple data for a fit
    :return x, y and e data
    """
    return (np.array([0, 1, 2, 3]),
            np.array([0.1, 1.2, 1.9, 3.15]),
            np.array([0.1, 0.09, 0.11, 0.1]))


def spline_data() -> (ndarray, ndarray, ndarray):
    """
    Gets 2 sets of data (high and low sampling)
    :return (high sampling:) x, y, e, (low sampling) x, y, e
    """
    x_data = np.linspace(0, 1, 20)
    np.random.seed(1)
    noise_stdev = 0.1
    y_data = np.random.normal(func(x_data), noise_stdev)
    e_data = 0.5*noise_stdev*np.ones(x_data.shape)
    # create low stats data
    xx = []
    yy = []
    ee = []

    for k in range(0, len(x_data), 2):
        xx.append(x_data[k])
        yy.append(y_data[k])
        ee.append(e_data[k])

    return x_data, y_data, e_data, np.array(xx), np.array(yy), np.array(ee)
