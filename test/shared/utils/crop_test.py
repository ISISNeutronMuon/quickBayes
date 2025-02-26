import unittest
import numpy as np
from quickBayes.utils.crop_data import crop


class CropTest(unittest.TestCase):

    def test_crop(self):
        x = np.linspace(0., 5., 11)
        y = x*2.
        e = x/10.

        x_crop, y_crop, e_crop = crop(x, y, e, 2.3, 4.8)

        expect_x = [2.5, 3.0, 3.5, 4.0, 4.5]
        expect_y = [5.0, 6.0, 7.0, 8.0, 9.0]
        expect_e = [.25, .30, .35, .40, .45]

        for j in range(len(x_crop)):
            self.assertAlmostEqual(x_crop[j], expect_x[j], 3)
            self.assertAlmostEqual(y_crop[j], expect_y[j], 3)
            self.assertAlmostEqual(e_crop[j], expect_e[j], 3)

    def test_crop_no_errors(self):
        x = np.linspace(0., 5., 11)
        y = x*2.
        e = None

        x_crop, y_crop, e_crop = crop(x, y, e, 2.3, 4.8)

        expect_x = [2.5, 3.0, 3.5, 4.0, 4.5]
        expect_y = [5.0, 6.0, 7.0, 8.0, 9.0]

        self.assertEqual(e_crop, None, 3)
        for j in range(len(x_crop)):
            self.assertAlmostEqual(x_crop[j], expect_x[j], 3)
            self.assertAlmostEqual(y_crop[j], expect_y[j], 3)

    def test_crop_start_before_data(self):
        x = np.linspace(0., 5., 11)
        y = x*2.
        e = x/10.

        x_crop, y_crop, e_crop = crop(x, y, e, -2.3, 4.8)

        expect_x = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        expect_y = [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        expect_e = [0, .05, 0.1, .15, 0.2, .25, .30, .35, .40, .45]

        for j in range(len(x_crop)):
            self.assertAlmostEqual(x_crop[j], expect_x[j], 3)
            self.assertAlmostEqual(y_crop[j], expect_y[j], 3)
            self.assertAlmostEqual(e_crop[j], expect_e[j], 3)

    def test_crop_end_after_data(self):
        x = np.linspace(0., 5., 11)
        y = x*2.
        e = x/10.

        x_crop, y_crop, e_crop = crop(x, y, e, 2.3, 7.8)

        expect_x = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        expect_y = [5.0, 6.0, 7.0, 8.0, 9.0, 10.]
        expect_e = [.25, .30, .35, .40, .45, 0.5]

        for j in range(len(x_crop)):
            self.assertAlmostEqual(x_crop[j], expect_x[j], 3)
            self.assertAlmostEqual(y_crop[j], expect_y[j], 3)
            self.assertAlmostEqual(e_crop[j], expect_e[j], 3)


if __name__ == '__main__':
    unittest.main()
