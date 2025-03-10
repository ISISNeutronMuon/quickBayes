import unittest
from quickBayes.functions.BG import FlatBG
from quickBayes.test_helpers.model_selection import WorkflowTemplateTest
from quickBayes.test_helpers.workflow_helper import gen_model_selection_data


class WorkflowTemplateWithGoFitTest(WorkflowTemplateTest,
                                    unittest.TestCase):

    def test_set_gofit_engine(self):
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_gofit_engine(5, [], [])
        self.assertEqual(self.wf.fit_engine.name, 'gofit')

    def test_update_gofit_engine(self):
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_gofit_engine(5, [], [])
        self.assertEqual(self.wf.fit_engine._lower, [])
        self.assertEqual(self.wf.fit_engine._upper, [])

        bg = FlatBG()
        self.wf.update_fit_engine(bg, [5])
        self.assertEqual(self.wf.fit_engine._lower, [-1])
        self.assertEqual(self.wf.fit_engine._upper, [1])


if __name__ == '__main__':
    unittest.main()
