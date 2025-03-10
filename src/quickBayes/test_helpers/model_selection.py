from quickBayes.workflow.model_selection.template import ModelSelectionWorkflow
from quickBayes.functions.BG import FlatBG
from quickBayes.functions.exp_decay import ExpDecay
from quickBayes.functions.composite import CompositeFunction
from quickBayes.test_helpers.workflow_helper import gen_model_selection_data


class SimpleWorkflow(ModelSelectionWorkflow):
    @staticmethod
    def _update_function(func):
        """
        Add a flat BG term
        """
        ed = ExpDecay()
        func.add_function(ed)
        return func


class WorkflowTemplateTest(object):

    def setUp(self):
        self.func = CompositeFunction()
        lin = FlatBG()
        self.func.add_function(lin)
        results = {}
        errors = {}
        self.wf = SimpleWorkflow(results, errors)

    def test_preprocess_data(self):
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.assertEqual(len(self.wf._data), 3)

        self.assertEqual(len(self.wf._raw['x']), len(x))
        self.assertEqual(len(self.wf._raw['y']), len(y))
        self.assertEqual(len(self.wf._raw['e']), len(e))

        self.assertEqual(len(self.wf._data['x']), len(x))
        self.assertEqual(len(self.wf._data['y']), len(y))
        self.assertEqual(len(self.wf._data['e']), len(e))
        for j in range(len(x)):
            self.assertEqual(self.wf._raw['x'][j], x[j])
            self.assertEqual(self.wf._raw['y'][j], y[j])
            self.assertEqual(self.wf._raw['e'][j], e[j])

            self.assertEqual(self.wf._data['x'][j], x[j])
            self.assertEqual(self.wf._data['y'][j], y[j])
            self.assertEqual(self.wf._data['e'][j], e[j])

    def test_get_raw_data(self):
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.assertEqual(len(self.wf._data), 3)

        raw = self.wf.get_raw
        self.assertEqual(len(raw['x']), len(x))
        self.assertEqual(len(raw['y']), len(y))
        self.assertEqual(len(raw['e']), len(e))

        for j in range(len(x)):
            self.assertEqual(raw['x'][j], x[j])
            self.assertEqual(raw['y'][j], y[j])
            self.assertEqual(raw['e'][j], e[j])

    def test_update_function(self):
        """
        Going to use get_guess to tell if the
        function has had anything added to it
        """
        self.assertEqual(len(self.func.get_guess()), 1)
        self.assertEqual(self.func._prefix, '')
        func = self.wf.update_function(self.func, 1)
        self.assertEqual(len(self.func.get_guess()), 3)
        for j, function in enumerate(func._funcs):
            self.assertEqual(function._prefix, f'N1:f{j + 1}.')

    def test_get_parameters_and_errors(self):
        params, errors = self.wf.get_parameters_and_errors
        self.assertEqual(params, {})
        self.assertEqual(errors, {})

        # setup workflow + generate data
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([0], [-9], [9])
        _ = self.wf.execute(1, self.func)
        # check parameters and errors updated
        params, errors = self.wf.get_parameters_and_errors
        self.assertEqual(len(params), 4)
        self.assertEqual(len(errors), 3)
        expected_keys = ['N1:f1.BG constant',
                         'N1:f2.Amplitude',
                         'N1:f2.lambda']

        self.assertEqual(list(errors.keys()), expected_keys)
        expected_param = [0.514, 1.000, 2.110]
        expected_error = [0.008, 0.913, 0.802]
        for j, key in enumerate(expected_keys):
            self.assertAlmostEqual(params[key][0], expected_param[j], 3)
            self.assertAlmostEqual(errors[key][0], expected_error[j], 3)

        prob_key = 'N1:loglikelihood'
        expected_keys.append(prob_key)
        self.assertEqual(list(params.keys()), expected_keys)
        self.assertAlmostEqual(params[prob_key][0], -72., 0)

    def test_fails_if_no_data(self):
        with self.assertRaises(ValueError):
            self.wf.set_scipy_engine([0], [-9], [9])

    def test_execute_no_engine(self):
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        with self.assertRaises(ValueError):
            _ = self.wf.execute(1, self.func)

    def test_add_second_engine_errors(self):
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([0], [-9], [9])
        with self.assertRaises(RuntimeError):
            self.wf.set_scipy_engine([0], [-9], [9])

    def test_set_scipy_engine(self):
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([], [], [])
        self.assertEqual(self.wf.fit_engine.name, 'scipy')

    def test_update_scipy_fit_engine(self):
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_model_selection_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([], [], [])
        self.assertEqual(self.wf.fit_engine._guess, [])
        self.assertEqual(self.wf.fit_engine._lower, [])
        self.assertEqual(self.wf.fit_engine._upper, [])

        bg = FlatBG()
        self.wf.update_fit_engine(bg, [2])
        self.assertEqual(self.wf.fit_engine._guess, [2])
        self.assertEqual(self.wf.fit_engine._lower, [-1])
        self.assertEqual(self.wf.fit_engine._upper, [1])
