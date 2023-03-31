from numpy import ndarray
import numpy as np
from abc import abstractmethod
from quickBayes.functions.BG import LinearBG
from quickBayes.test_helpers.fitting_data import (basic_data,
                                                  spline_data,
                                                  func)


"""
This is a testing template for the fit engines.
This provides the minimum amount of testing.
If a fit engine has extra features they will
need to be tested too.
"""


class FitEngineTemplate(object):
    """
    This class gives a basic outline for testing
    a fitting engine.

    ---------IMPORTANT-------------------------------
    When inheriting this class, you will also need to
    inherit "unittest.TestCase" for it to work as a unit
    test. Cannot do it here as it will then treat this base
    class as a test (which will always fail).
    -------------------------------------------------

    The below methods give the expected results
    These are the only bits that need changing for a
    new fit engine instance.
    """
    @abstractmethod
    def get_test_engine(self, x: ndarray, y: ndarray, e: ndarray):
        raise NotImplementedError()

    @abstractmethod
    def get_name(self) -> str:
        """
        :return the expected name
        """
        raise NotImplementedError()

    @abstractmethod
    def get_basic_fit_params(self) -> (ndarray, ndarray):
        """
        :return the expected parameters and errors
        """
        raise NotImplementedError()

    @abstractmethod
    def get_basic_fit_values(self) -> (ndarray, ndarray,
                                       ndarray, ndarray):
        """
        :return the expected; y, error, difference and difference errors
        """
        raise NotImplementedError()

    @abstractmethod
    def get_chi_squared(self) -> float:
        """
        :return the expected chi^2 value
        """
        raise NotImplementedError()

    @abstractmethod
    def get_covariance(self) -> ndarray:
        """
        :return the expected covariance matrix
        """
        raise NotImplementedError()

    @abstractmethod
    def get_spline_params(self) -> (ndarray, ndarray):
        """
        :return the expected parameters and errors for splined data
        """
        raise NotImplementedError()

    @abstractmethod
    def get_spline_fits(self) -> (ndarray, ndarray,
                                  ndarray, ndarray):
        """
        :return the expected; y, e, differences, difference errors
        for the splined data
        """
        raise NotImplementedError()

    @abstractmethod
    def get_low_stat_params(self) -> (ndarray, ndarray):
        """
        :return the expected parameters and errors
        """
        raise NotImplementedError()

    @abstractmethod
    def get_low_stat_fits(self) -> (ndarray, ndarray,
                                    ndarray, ndarray):
        """
        :return the expected; y, e, differences, difference errors
        """
        raise NotImplementedError()

    @abstractmethod
    def get_spline_chi2(self) -> dict:
        """
        :return a dict of the chi^2 values (keys= low, high)
        """
        raise NotImplementedError()

    @abstractmethod
    def get_spline_covar(self) -> dict:
        """
        :return a dict of the covar matricies (keys= low, high)
        """
        raise NotImplementedError()

    """
    Below are the helpers for asserting data
    """
    def assert_parameters(self, params: ndarray, errors: ndarray,
                          expected_p: ndarray,
                          expected_e: ndarray) -> None:
        """
        Method for asserting that the calculated and expected
        parameters match
        :param params: the parameters calculated from fit
        :param errors: the errors calculated from fit
        :param expected_p: the expected parameters
        :param expected_e: the expected parameter errors
        """
        self.assertEqual(len(params), len(errors))
        for k in range(len(params)):
            self.assertAlmostEqual(params[k], expected_p[k], 3)
            self.assertAlmostEqual(errors[k], expected_e[k], 3)

    def assert_fit_values(self, xf: ndarray, yf: ndarray, ef: ndarray,
                          df: ndarray, de: ndarray, x_data: ndarray,
                          expected_y: ndarray, expected_e: ndarray,
                          expected_d: ndarray,
                          expected_de: ndarray) -> None:
        """
        Method for testing that the fit values and expected
        results match.
        :param xf: the x data from the fit
        :param yf: the y data from the fit
        :param ef: the e data from the fit
        :param df: the difference data from the fit
        :param de: the difference errors from the fit
        :param x_data: the x data that should have been used
        :param expected_y: the expected y data
        :param expected_e: the expected e data
        :param expected_d: the expected difference data
        :param expected_de: the expected difference errors

        """

        self.assertEqual(len(xf), len(x_data))
        self.assertEqual(len(yf), len(x_data))
        self.assertEqual(len(ef), len(x_data))
        self.assertEqual(len(df), len(x_data))
        self.assertEqual(len(de), len(x_data))
        for k in range(len(xf)):
            self.assertAlmostEqual(xf[k], x_data[k], 3)
            self.assertAlmostEqual(yf[k], expected_y[k], 3)
            self.assertAlmostEqual(ef[k], expected_e[k], 3)
            self.assertAlmostEqual(df[k], expected_d[k], 3)
            self.assertAlmostEqual(de[k], expected_de[k], 3)

    def assert_covar_matrix(self, covar: ndarray,
                            expect: ndarray) -> None:
        """
        Method to test that 2D matrices match
        :param covar: calculated covariance matrix
        :param expect: the expected covariance matrix
        """
        # assume a 2D matrix
        self.assertEqual(len(covar), len(expect))
        self.assertEqual(len(covar[0]), len(expect[0]))

        for i in range(len(covar)):
            for j in range(len(covar[i])):
                self.assertAlmostEqual(covar[i][j], expect[i][j], 3)

    """
    Below are the default fitting tests
    """
    def test_name(self) -> None:
        """
        Test the fit engines name
        """
        x_data, y_data, e_data = basic_data()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        expect = self.get_name()
        self.assertEqual(self.engine.name, expect)

    def test_fit_params(self) -> None:
        """
        Test the fit engine gets the correct fit parameters
        """
        x_data, y_data, e_data = basic_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)

        params, errors = self.engine.get_fit_parameters()
        expected_p, expected_e = self.get_basic_fit_params()

        self.assert_parameters(params, errors, expected_p,
                               expected_e)

    def test_fit_values(self) -> None:
        """
        Test the fit engine gets the expected fit values
        """
        x_data, y_data, e_data = basic_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)

        xf, yf, ef, df, de = self.engine.get_fit_values()
        (expected_y, expected_e, expected_d,
         expected_de) = self.get_basic_fit_values()
        self.assert_fit_values(xf, yf, ef, df, de, x_data,
                               expected_y, expected_e,
                               expected_d, expected_de)

    def test_chi_squared(self) -> None:
        """
        Test the fit engine gets the expected chi^2
        """
        x_data, y_data, e_data = basic_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)
        expect = self.get_chi_squared()
        self.assertAlmostEqual(self.engine.get_chi_squared(), expect, 3)

    def test_cov(self) -> None:
        """
        Test that the fit engine gets the expected covariance matrix
        Need to do more data to get an accurate covariance matrix
        """
        x_data = np.linspace(0, 1, 10)
        y_data = func(x_data)
        e_data = 0.1*np.ones(len(x_data))
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)
        calculated = self.engine.get_covariance_matrix()
        expect = self.get_covariance()

        self.assert_covar_matrix(calculated, expect)

    def fit_data_with_diff_sampling(self):
        x_data, y_data, e_data, xx, yy, ee = spline_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(xx, yy, ee)

        # fit with less data
        self.engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        self.engine.do_fit(x_data, y_data, e_data, bg)
        return xx

    def test_spline_data_params(self) -> None:
        """
        Want to test fitting data with high and low sampling/stats.
        These should both be recorded (history).
        Both should match the original data sampling.
        This test is for the parameters.
        """
        _ = self.fit_data_with_diff_sampling()
        # check latest results
        params, errors = self.engine.get_fit_parameters()
        expected_p, expected_e = self.get_spline_params()
        self.assert_parameters(params, errors, expected_p, expected_e)

        # check first results -> lower stats
        params, errors = self.engine.get_fit_parameters(0)
        expected_p, expected_e = self.get_low_stat_params()
        self.assert_parameters(params, errors, expected_p, expected_e)

    def test_spline_data_fits(self) -> None:
        """
        Want to test fitting data with high and low sampling/stats.
        These should both be recorded (history).
        Both should match the original data sampling.
        This test is for the fits.
        """
        xx = self.fit_data_with_diff_sampling()

        # check latest results
        xf, yf, ef, df, de = self.engine.get_fit_values()
        (expected_y, expected_e,
         expected_d, expected_de) = self.get_spline_fits()

        self.assert_fit_values(xf, yf, ef, df, de, xx, expected_y,
                               expected_e, expected_d, expected_de)

        # check first results -> lower stats
        xf, yf, ef, df, de = self.engine.get_fit_values(0)
        (expected_y, expected_e,
         expected_d, expected_de) = self.get_low_stat_fits()
        self.assert_fit_values(xf, yf, ef, df, de, xx, expected_y,
                               expected_e, expected_d, expected_de)

    def test_spline_chi_squared(self):
        """
        Want to test fitting data with high and low sampling/stats.
        These should both be recorded (history).
        Both should match the original data sampling.
        This test is for the chi^2.
        """
        _ = self.fit_data_with_diff_sampling()

        expect = self.get_spline_chi2()
        self.assertAlmostEqual(self.engine.get_chi_squared(),
                               expect['high'], 3)
        self.assertAlmostEqual(self.engine.get_chi_squared(0),
                               expect['low'], 3)

    def test_spline_cov(self):
        """
        Want to test fitting data with high and low sampling/stats.
        These should both be recorded (history).
        Both should match the original data sampling.
        This test is for the covariance matrices.
        """
        _ = self.fit_data_with_diff_sampling()

        calculated = self.engine.get_covariance_matrix()
        expect = self.get_spline_covar()

        self.assert_covar_matrix(calculated, expect['high'])

        # check 1st fit
        calculated = self.engine.get_covariance_matrix(0)
        self.assert_covar_matrix(calculated, expect['low'])
