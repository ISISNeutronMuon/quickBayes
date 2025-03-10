import numpy as np
from quickBayes.test_helpers.template_fit_test import FitEngineTemplate


class ScipyFitTemplate(FitEngineTemplate):

    @staticmethod
    def get_test_engine(x, y, e):
        raise NotImplementedError()

    @staticmethod
    def get_name():
        raise NotImplementedError()

    @staticmethod
    def get_basic_fit_params():
        raise NotImplementedError()

    @staticmethod
    def get_chi_squared():
        return 2.347

    @staticmethod
    def get_covariance():
        raise NotImplementedError()

    @staticmethod
    def get_basic_fit_values():
        raise NotImplementedError()

    @staticmethod
    def get_spline_params():
        return [0.882, 0.095], [0.037, 0.022]

    @staticmethod
    def get_spline_fits():
        expected_y = [.095, 0.188, 0.281, 0.374, 0.467, 0.559,
                      0.652, 0.745, 0.838, 0.931]
        expected_e = [0.022, 0.019, 0.016, 0.013, 0.012, 0.012,
                      0.013, 0.015, 0.017, 0.020]
        expected_d = [-0.167, 0.046, -0.095, -0.185, -0.044, -0.161,
                      0.016, -0.132, -0.003, -0.026]
        expected_de = [0.055, 0.053, 0.052, 0.052, 0.051, 0.051,
                       0.052, 0.052, 0.053, 0.054]
        return expected_y, expected_e, expected_d, expected_de

    @staticmethod
    def get_low_stat_params():
        return [0.811, 0.204], [0.052, 0.029]

    @staticmethod
    def get_low_stat_fits():
        expected_y = [.204, 0.289, 0.375, 0.460, 0.545, 0.631,
                      0.716, 0.801, 0.887, 0.972]
        expected_e = [0.031, 0.027, 0.022, 0.019, 0.017, 0.017,
                      0.019, 0.022, 0.027, 0.031]
        expected_d = [-0.059, 0.147, -0.001, -0.099, 0.034, -0.089,
                      0.080, -0.075, 0.046, 0.015]
        expected_de = [0.059, 0.057, 0.055, 0.053, 0.053, 0.053,
                       0.053, 0.055, 0.057, 0.059]

        return expected_y, expected_e, expected_d, expected_de

    @staticmethod
    def get_spline_chi2():
        return {'low': 2.921, 'high': 5.332}

    @staticmethod
    def get_spline_covar():
        high = np.array([np.array([0.001, -0.001]),
                         np.array([-0.001, 0.0004])])
        low = np.array([np.array([0.003, -0.001]),
                        np.array([-0.001, 0.0004])])

        return {'high': high, 'low': low}
