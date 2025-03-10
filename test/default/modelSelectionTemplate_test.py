import unittest
from quickBayes.test_helpers.model_selection import WorkflowTemplateTest
from quickBayes.test_helpers.workflow_helper import gen_model_selection_data


class WorkflowTemplateWithoutGoFitTest(WorkflowTemplateTest,
                                       unittest.TestCase):

    def test_set_gofit_engine(self):
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        with self.assertRaises(RuntimeError):
            self.wf.set_gofit_engine(5, [], [])


if __name__ == '__main__':
    unittest.main()
