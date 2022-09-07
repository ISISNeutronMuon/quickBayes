"""Characterization tests for QLres module"""
import unittest
from quasielasticbayes.data import DatCom
from quasielasticbayes.util import LUDCMP  # noqa: F401


class StuffTest(unittest.TestCase):
    def test_primes(self):
        a = DatCom(2, 3)
        self.assertTrue(a is not None)
        self.assertEqual(1, 1)


if __name__ == '__main__':
    unittest.main()
