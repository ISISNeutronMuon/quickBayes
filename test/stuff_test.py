"""Characterization tests for QLres module"""
import unittest
from quasielasticbayes.data import DatCom


class StuffTest(unittest.TestCase):
    def test_primes(self):
        a = DatCom
        self.assert(a)
        self.assertEqual(1, 1)


if __name__ == '__main__':
    unittest.main()
