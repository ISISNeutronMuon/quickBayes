"""Characterization tests for QLres module"""
import unittest
from quasielasticbayes.stuff import primes


class StuffTest(unittest.TestCase):
    def test_primes(self):
        self.assertEqual([2, 3, 5], primes(3))


if __name__ == '__main__':
    unittest.main()
