import unittest
from quickBayes.utils.general import get_background_function


class GetBGTest(unittest.TestCase):
    """
    Since the 3 BG's all have different
    number of parameters, we can use
    that information to test if
    we have got the correct function.
    """

    def test_none(self):
        func = get_background_function("none")
        self.assertEqual(func.N_params, 0)

    def test_none_mixed_case(self):
        func = get_background_function("NonE")
        self.assertEqual(func.N_params, 0)

    def test_flat(self):
        func = get_background_function("flat")
        self.assertEqual(func.N_params, 1)

    def test_flat_mixed_case(self):
        func = get_background_function("FlAT")
        self.assertEqual(func.N_params, 1)

    def test_linear(self):
        func = get_background_function("linear")
        self.assertEqual(func.N_params, 2)

    def test_linear_mixed_case(self):
        func = get_background_function("LinEAR")
        self.assertEqual(func.N_params, 2)

    def test_invalid(self):
        with self.assertRaises(ValueError):
            self.assertRaises(get_background_function("rubbish"))


if __name__ == '__main__':
    unittest.main()
