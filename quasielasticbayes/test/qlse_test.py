"""Characterization tests for QLres module"""
import os.path
import unittest
import numpy as np
from quasielasticbayes.testing import load_json, add_path
import tempfile
from quasielasticbayes.QLres import qlres
import json
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
class QLresTest(unittest.TestCase):
    """
    Characterization tests using inputs that have been accepted as correct.
    The output is based on running BayesQuasi algorithm in mantid 6.2
    with the inputs taken from the BayesQuasiTest unit test
    """

    def test_qlres_minimal_input(self):
        # reference inputs
        with open(os.path.join(DATA_DIR, 'qlse', 'qlse_input.json'), 'r') as fh:
            inputs = load_json(fh)
        temp_dir = tempfile.TemporaryDirectory()	
        inputs['wrks'] = add_path(temp_dir.name, inputs['wrks'])
        nd, xout, yout, eout, yfit, yprob = qlres(inputs['numb'], inputs['Xv'], inputs['Yv'], inputs['Ev'],
                                                  inputs['reals'], inputs['fitOp'],
                                                  inputs['Xdat'], inputs['Xb'], inputs['Yb'],
                                                  inputs['Wy'], inputs['We'], inputs['dtn'], inputs['xsc'],
                                                  inputs['wrks'], inputs['wrkr'], inputs['lwrk'])

        # verify
        with open(os.path.join(DATA_DIR, 'qlse', 'qlse_output.json'), 'r') as fh:
            reference = load_json(fh)
        
        self.assertEqual(reference['nd'], nd)
        np.testing.assert_almost_equal(reference['xout'], xout)
        np.testing.assert_almost_equal(reference['yout'], yout)
        np.testing.assert_almost_equal(reference['eout'], eout)
        np.testing.assert_almost_equal(reference['yfit'], yfit)
        np.testing.assert_almost_equal(reference['yprob'], yprob)
        temp_dir.cleanup()
if __name__ == '__main__':
    unittest.main()