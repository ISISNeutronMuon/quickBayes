"""Characterization tests for QLres module"""
import os.path
import unittest
import numpy as np
from quasielasticbayes.testing import load_json, add_path

from quasielasticbayes.QLres import qlres

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class QLresTest(unittest.TestCase):
    """
    Characterization tests using inputs that have been accepted as correct.
    The output is based on running BayesQuasi algorithm in mantid 6.2
    with the inputs taken from the BayesQuasiTest unit test
    """

    def test_qlres_minimal_input(self):
        # reference inputs
        with open(os.path.join(DATA_DIR, 'qlres', 'qlres-input-spec-0.json'), 'r') as fh:
            inputs = load_json(fh)
        inputs['wrks'] = add_path(inputs['wrks'])
        nd, xout, yout, eout, yfit, yprob = qlres(inputs['numb'], inputs['Xv'], inputs['Yv'], inputs['Ev'],
                                                  inputs['reals'], inputs['fitOp'],
                                                  inputs['Xdat'], inputs['Xb'], inputs['Yb'],
                                                  inputs['Wy'], inputs['We'], inputs['dtn'], inputs['xsc'],
                                                  inputs['wrks'], inputs['wrkr'], inputs['lwrk'])

        # verify
        with open(os.path.join(DATA_DIR, 'qlres', 'qlres-output-spec-0.json'), 'r') as fh:
            reference = load_json(fh)
        
        self.assertEqual(reference['nd'], nd)
        np.testing.assert_almost_equal(reference['xout'], xout)
        np.testing.assert_almost_equal(reference['yout'], yout)
        np.testing.assert_almost_equal(reference['eout'], eout)
        np.testing.assert_almost_equal(reference['yfit'], yfit)
        np.testing.assert_almost_equal(reference['yprob'], yprob)

if __name__ == '__main__':
    unittest.main()