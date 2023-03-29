from quasielasticbayes.v2.workflow.MuonExpDecay import muon_expdecay_main
from quasielasticbayes.v2.workflow.QlData import ql_data_main
from quasielasticbayes.v2.functions.qse_fixed import QSEFixFunction
from quasielasticbayes.v2.functions.BG import LinearBG
from quasielasticbayes.v2.workflow.qse_search import QSEGridSearch
from quasielasticbayes.v2.workflow.QSE import qse_data_main
import numpy as np
import time
import os.path
import unittest

DATA_DIR = os.path.join(os.path.dirname(__file__),  'data')
MUON_DIR = os.path.join(DATA_DIR, 'muon')


def run_muon():
    sx, sy, se = np.loadtxt(os.path.join(MUON_DIR,
                                             'muon_expdecay_3_big.npy'))
    sample = {'x': sx, 'y': sy, 'e': se}
    results = {}
    errors = {}
    start = time.time()
    (results, errors,
     new_x, fits, f_errors) = muon_expdecay_main(sample, "flat",
                                                 0.11, 15.0,
                                                 results, errors)
    return time.time() - start


def run_qldata():
    sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
    rx, ry, re = np.load(os.path.join(DATA_DIR, 'resolution_data_red.npy'))

    sample = {'x': sx, 'y': sy, 'e': se}
    resolution = {'x': rx, 'y': ry}
    results = {}
    errors = {}
    start = time.time()
    (results, errors,
     new_x, fits, f_errors) = ql_data_main(sample, resolution,
                                           "linear", -0.4, 0.4,
                                           True, results, errors)
    return time.time() - start


def run_qse():
    sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
    rx, ry, re = np.load(os.path.join(DATA_DIR, 'qse_res.npy'),
                         allow_pickle=True)

    sample = {'x': sx, 'y': sy, 'e': se}
    resolution = {'x': rx, 'y': ry}
    results = {}
    errors = {}
    start = time.time()
    (results, errors,
     new_x, fits, fit_e) = qse_data_main(sample, resolution,
                                         "linear", -0.4, 0.4, True,
                                         results, errors)
    return time.time() - start


def run_quest():
    sx, sy, se = np.load(os.path.join(DATA_DIR, 'sample_data_red.npy'))
    rx, ry, re = np.load(os.path.join(DATA_DIR, 'qse_res.npy'),
                             allow_pickle=True)
    resolution = {'x': rx, 'y': ry}

    # set values
    start_x = -0.5
    end_x = 0.5

    N_beta = 5
    beta_start = 0.75
    beta_end = 0.78

    N_FWHM = 5
    FWHM_start = 0.055
    FWHM_end = 0.056

    # create workflow
    start = time.time()
    search = QSEGridSearch()
    new_x, ry = search.preprocess_data(sx, sy, se, start_x, end_x,
                                       resolution)
    search.set_x_axis(beta_start, beta_end, N_beta, 'beta')
    search.set_y_axis(FWHM_start, FWHM_end, N_FWHM, 'FWHM')

    # create function
    bg = LinearBG()
    func = QSEFixFunction(bg, True, new_x, ry, start_x, end_x)
    func.add_single_SE()
    func.set_delta_bounds([0, -.5], [20, .5])

    # do search
    search.set_scipy_engine(func.get_guess(), *func.get_bounds())
    X, Y = search.execute(func)
    return time.time() - start


class proTest(unittest.TestCase):

    def test_muon(self):
        times = []
        for j in range(10):
            times.append(run_muon())
        times = np.array(times)
        print("\n muon", np.mean(times), np.std(times))
        # self.assertEqual(1, 2)

    def test_qldata(self):
        times = []
        for j in range(10):
            times.append(run_qldata())
        times = np.array(times)
        print("\n qldata", np.mean(times), np.std(times))

        # self.assertEqual(1, 2)

    def test_qlse(self):
        times = []
        for j in range(10):
            times.append(run_qse())
        times = np.array(times)
        print("\n qse", np.mean(times), np.std(times))
        # self.assertEqual(1, 2)

    def test_quest(self):
        times = []
        for j in range(10):
            times.append(run_quest())
        times = np.array(times)
        print("\n quest", np.mean(times), np.std(times))
        # self.assertEqual(1, 2)


if __name__ == '__main__':
    unittest.main()
