from numpy import ndarray
import numpy as np
from abc import abstractmethod
from quasielasticbayes.v2.functions.BG import LinearBG
from quasielasticbayes.test_helpers.fitting_data import (basic_data,
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
    The below methods give the expected results
    These are the only bits that need changing for a
    new fit engine instance.
    """
    @abstractmethod
    def get_test_engine(self, x: ndarray, y: ndarray, e: ndarray):
        raise NotImplementedError()

    @abstractmethod
    def get_name(self):
        raise NotImplementedError()

    @abstractmethod
    def get_basic_fit_params(self):
        raise NotImplementedError()

    @abstractmethod
    def get_basic_fit_values(self):
        raise NotImplementedError()

    @abstractmethod
    def get_chi_squared(self):
        raise NotImplementedError()

    @abstractmethod
    def get_covariance(self):
        raise NotImplementedError()

    @abstractmethod
    def get_spline_params(self):
        raise NotImplementedError()

    @abstractmethod
    def get_spline_fits(self):
        raise NotImplementedError()

    @abstractmethod
    def get_low_stat_params(self):
        raise NotImplementedError()

    @abstractmethod
    def get_low_stat_fits(self):
        raise NotImplementedError()

    @abstractmethod
    def get_spline_chi2(self):
        raise NotImplementedError()

    @abstractmethod
    def get_spline_covar(self):
        raise NotImplementedError()

    """
    Below are the helpers for asserting data
    """
    def assert_parameters(self, params, errors, expected_p,
                          expected_e):
        self.assertEqual(len(params), len(errors))
        for k in range(len(params)):
            self.assertAlmostEqual(params[k], expected_p[k], 3)
            self.assertAlmostEqual(errors[k], expected_e[k], 3)

    def assert_fit_values(self, xf, yf, ef, df, de, x_data,
                          expected_y, expected_e, expected_d, expected_de):

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

    def assert_covar_matrix(self, covar, expect):
        # assume a 2D matrix
        self.assertEqual(len(covar), len(expect))
        self.assertEqual(len(covar[0]), len(expect[0]))

        for i in range(len(covar)):
            for j in range(len(covar[i])):
                self.assertAlmostEqual(covar[i][j], expect[i][j], 3)

    """
    Below are the default fitting tests
    """
    def test_name(self):
        x_data, y_data, e_data = basic_data()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        expect = self.get_name()
        self.assertEqual(self.engine.name, expect)

    def test_fit_params(self):
        x_data, y_data, e_data = basic_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)

        params, errors = self.engine.get_fit_parameters()
        expected_p, expected_e = self.get_basic_fit_params()

        self.assert_parameters(params, errors, expected_p,
                               expected_e)

    def test_fit_values(self):
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

    def test_chi_squared(self):
        x_data, y_data, e_data = basic_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)
        expect = self.get_chi_squared()
        self.assertAlmostEqual(self.engine.get_chi_squared(), expect, 3)

    def test_cov(self):
        # need to do more data to get an accurate covariance matrix
        x_data = np.linspace(0, 1, 10)
        y_data = func(x_data)
        e_data = 0.1*np.ones(len(x_data))
        bg = LinearBG()

        self.engine = self.get_test_engine(x_data, y_data, e_data)
        self.engine.do_fit(x_data, y_data, e_data, bg)
        calculated = self.engine.get_covariance_matrix()
        expect = self.get_covariance()

        self.assert_covar_matrix(calculated, expect)

    def test_spline_data_params(self):
        x_data, y_data, e_data, xx, yy, ee = spline_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(xx, yy, ee)

        # fit with less data
        self.engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        self.engine.do_fit(x_data, y_data, e_data, bg)

        # check latest results
        params, errors = self.engine.get_fit_parameters()
        expected_p, expected_e = self.get_spline_params()

        self.assert_parameters(params, errors, expected_p, expected_e)

        # check first results -> lower stats
        params, errors = self.engine.get_fit_parameters(0)
        expected_p, expected_e = self.get_low_stat_params()

        self.assert_parameters(params, errors, expected_p, expected_e)

    def test_spline_data_fits(self):
        x_data, y_data, e_data, xx, yy, ee = spline_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(xx, yy, ee)

        # fit with less data
        self.engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        self.engine.do_fit(x_data, y_data, e_data, bg)

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
        x_data, y_data, e_data, xx, yy, ee = spline_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(xx, yy, ee)

        # fit with less data
        self.engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        self.engine.do_fit(x_data, y_data, e_data, bg)

        expect = self.get_spline_chi2()
        self.assertAlmostEqual(self.engine.get_chi_squared(),
                               expect['high'], 3)
        self.assertAlmostEqual(self.engine.get_chi_squared(0),
                               expect['low'], 3)

    def test_spline_cov(self):
        x_data, y_data, e_data, xx, yy, ee = spline_data()
        bg = LinearBG()

        self.engine = self.get_test_engine(xx, yy, ee)

        # fit with less data
        self.engine.do_fit(xx, yy, ee, bg)
        # fit with more data
        self.engine.do_fit(x_data, y_data, e_data, bg)

        calculated = self.engine.get_covariance_matrix()
        expect = self.get_spline_covar()

        self.assert_covar_matrix(calculated, expect['high'])

        # check 1st fit
        calculated = self.engine.get_covariance_matrix(0)
        self.assert_covar_matrix(calculated, expect['low'])
