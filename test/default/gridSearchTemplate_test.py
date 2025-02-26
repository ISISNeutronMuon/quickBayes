import unittest
from quickBayes.test_helpers.workflow_helper import gen_grid_search_data
from quickBayes.test_helpers.grid_template import GridSearchTemplateTest


class GridSearchTemplateWithoutGoFitTest(GridSearchTemplateTest,
                                         unittest.TestCase):
    def test_set_gofit_engine(self):
        # same
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        with self.assertRaises(RuntimeError):
            self.wf.set_gofit_engine(5, [], [])


if __name__ == '__main__':
    unittest.main()
