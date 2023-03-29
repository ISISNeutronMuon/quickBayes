import unittest
from quickBayes.utils.general import update_guess
from quickBayes.functions.BG import LinearBG
from quickBayes.functions.BG import NoBG


class UpdateGuessTest(unittest.TestCase):

    """
    This expects the whole function.
    So it will not update params from
    the function it already has
    """

    def test_update_guess(self):
        params = [1]
        func = LinearBG()
        result = update_guess(params, func)
        expect = [1, 0]
        self.assertEqual(len(result), len(expect))
        for j in range(len(expect)):
            self.assertEqual(result[j], expect[j])

    def test_update_guess_no_new_guess(self):
        params = [1, 2]
        func = LinearBG()
        result = update_guess(params, func)
        self.assertEqual(len(result), len(params))
        for j in range(len(params)):
            self.assertEqual(result[j], params[j])

    def test_update_guess_all_empty(self):
        params = []
        func = NoBG()
        result = update_guess(params, func)
        expect = []
        self.assertEqual(len(result), len(expect))

    def test_update_guess_empty_params(self):
        params = []
        func = LinearBG()
        result = update_guess(params, func)
        expect = [0, 0]
        self.assertEqual(len(result), len(expect))
        for j in range(len(expect)):
            self.assertEqual(result[j], expect[j])

    def test_invalid(self):
        params = [1, 2, 3]
        func = LinearBG()
        with self.assertRaises(ValueError):
            _ = update_guess(params, func)


if __name__ == '__main__':
    unittest.main()
