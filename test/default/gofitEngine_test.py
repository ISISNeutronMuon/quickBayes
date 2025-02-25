import unittest
from quickBayes.fitting.gofit_engine import GoFitEngine


class GoFitEngineNotInstalledTest(unittest.TestCase):
    def test_gofit(self):
        with self.assertRaises(RuntimeError):
            GoFitEngine([1], [1], [1],
                        [1], [1], 1, 1)


if __name__ == '__main__':
    unittest.main()
