"""Characterization tests for QLres module"""
import os.path
import unittest
import numpy as np
from quasielasticbayes.testing import load_json, add_path
import tempfile
from quasielasticbayes.QLdata import qldata
from quasielasticbayes.python.qldata_main import QLdata

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class QLdataTest(unittest.TestCase):
    """
    Characterization tests using inputs that have been accepted as correct.
    The output is based on running BayesQuasi algorithm in mantid 6.2
    with the inputs taken from the BayesQuasiTest unit test
    """

    def test_qlres_minimal_input(self):
        # reference inputs
        with open(os.path.join(DATA_DIR, 'qldata', 'qldata_input.json'), 'r') as fh:
            inputs = load_json(fh)
        temp_dir = os.path.join(os.path.realpath(__file__), "..","..") #tempfile.TemporaryDirectory()	
        inputs['wrks'] = add_path(temp_dir, inputs['wrks'])
        #inputs['wrks'] = add_path(temp_dir.name, inputs['wrks'])
        nd, xout, yout, eout, yfit, yprob = qldata(inputs['numb'], inputs['Xv'], inputs['Yv'], inputs['Ev'],
                                                  inputs['reals'], inputs['fitOp'],
                                                  inputs['Xdat'], inputs['Xb'], inputs['Yb'], inputs['Eb'],
                                                  inputs['Wy'], inputs['We'], 
                                                  inputs['wrks'], inputs['wrkr'], inputs['lwrk'])
        nd, xout, yout, eout, yfit, yprob = QLdata(inputs['numb'], inputs['Xv'], inputs['Yv'], inputs['Ev'],
                                                  inputs['reals'], inputs['fitOp'],
                                                  inputs['Xdat'], inputs['Xb'], inputs['Yb'], inputs['Eb'],
                                                  inputs['Wy'], inputs['We'], 
                                                  inputs['wrks'], inputs['wrkr'], inputs['lwrk'])
        # verify
        with open(os.path.join(DATA_DIR, 'qldata', 'qldata_output.json'), 'r') as fh:
            reference = load_json(fh)
        
        self.assertEqual(reference['nd'], nd)
        np.testing.assert_almost_equal(reference['xout'], xout)
        np.testing.assert_almost_equal(reference['yout'], yout)
        np.testing.assert_almost_equal(reference['eout'], eout)
        np.testing.assert_almost_equal(reference['yfit'], yfit)
        np.testing.assert_almost_equal(reference['yprob'], yprob)
        #temp_dir.cleanup()


if __name__ == '__main__':
    unittest.main()
