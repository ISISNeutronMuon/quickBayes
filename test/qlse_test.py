"""Characterization tests for QLres module"""
import os.path
import unittest
import numpy as np
from quasielasticbayes.testing import load_json, add_path
from quasielasticbayes.testing import get_OS_precision, get_qlse_prob
import tempfile
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
        fin = 'qlse_input.json'
        with open(os.path.join(DATA_DIR, 'qlse', fin), 'r') as fh:
            inputs = load_json(fh)
        with tempfile.TemporaryDirectory() as tmp_dir:

            inputs['wrks'] = add_path(tmp_dir, inputs['wrks'])
            nd, xout, yout, eout, yfit, yprob = qlres(inputs['numb'],
                                                      inputs['Xv'],
                                                      inputs['Yv'],
                                                      inputs['Ev'],
                                                      inputs['reals'],
                                                      inputs['fitOp'],
                                                      inputs['Xdat'],
                                                      inputs['Xb'],
                                                      inputs['Yb'],
                                                      inputs['Wy'],
                                                      inputs['We'],
                                                      inputs['dtn'],
                                                      inputs['xsc'],
                                                      inputs['wrks'],
                                                      inputs['wrkr'],
                                                      inputs['lwrk'])

            # verify
            cf = 'qlse_output.json'
            with open(os.path.join(DATA_DIR, 'qlse', cf), 'r') as fh:
                reference = load_json(fh)

            dp = get_OS_precision()
            self.assertEqual(reference['nd'], nd)
            np.testing.assert_almost_equal(reference['xout'], xout,
                                           decimal=dp)
            np.testing.assert_almost_equal(reference['yout'], yout,
                                           decimal=dp)
            np.testing.assert_almost_equal(reference['eout'], eout,
                                           decimal=dp)
            np.testing.assert_almost_equal(reference['yfit'], yfit,
                                           decimal=dp)
            ref_prob = get_qlse_prob(reference['yprob'])
            np.testing.assert_almost_equal(ref_prob, yprob,
                                           decimal=dp)


if __name__ == '__main__':
    unittest.main()
