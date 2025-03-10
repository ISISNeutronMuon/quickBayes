from quickBayes.workflow.grid_search.template import GridSearchTemplate
from quickBayes.functions.BG import FlatBG
from quickBayes.functions.exp_decay import ExpDecay
from quickBayes.test_helpers.workflow_helper import (gen_grid_search_data,
                                                     FixedBG,
                                                     FixedComposite)


class SimpleWorkflow(GridSearchTemplate):
    @staticmethod
    def _update_function(func):
        """
        Add a flat BG term
        """
        ed = ExpDecay()
        func.add_function(ed)
        return func

    @staticmethod
    def _set_x_value(func, value):
        func.set_c(value)
        return func

    @staticmethod
    def _set_y_value(func, value):
        func.set_m(value)
        return func

    @staticmethod
    def N(func):
        return 1


class GridSearchTemplateTest(object):

    def setUp(self):
        self.func = FixedComposite()
        lin = FixedBG()
        self.func.add_function(lin)
        self.wf = SimpleWorkflow()

    def test_set_x_axis(self):
        self.wf.set_x_axis(1, 3, 5, 'test')
        x_axis = self.wf.get_x_axis
        expect = [1, 1.5, 2., 2.5, 3]
        values = x_axis.values
        self.assertEqual(len(values), len(expect))
        for j in range(len(expect)):
            self.assertEqual(values[j], expect[j])
        self.assertEqual(x_axis.len, 5)
        self.assertEqual(x_axis.label, "test")

    def test_set_y_axis(self):
        self.wf.set_y_axis(6, 8, 5, 'test2')
        y_axis = self.wf.get_y_axis
        expect = [6, 6.5, 7., 7.5, 8]
        values = y_axis.values
        self.assertEqual(len(values), len(expect))
        for j in range(len(expect)):
            self.assertEqual(values[j], expect[j])
        self.assertEqual(y_axis.len, 5)
        self.assertEqual(y_axis.label, "test2")

    def test_grid(self):
        self.wf.set_x_axis(0, 1, 2, 'x')
        self.wf.set_y_axis(1, 2, 2, 'y')
        X, Y = self.wf._generate_grid()
        grid = self.wf.get_grid

        expect_x = [[0, 1], [0, 1]]
        expect_y = [[1, 1], [2, 2]]

        for i in range(2):
            for j in range(2):
                self.assertEqual(X[i][j], expect_x[i][j])
                self.assertEqual(Y[i][j], expect_y[i][j])
                self.assertEqual(grid[i][j], 0)

    def test_preprocess_data(self):
        # same
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.assertEqual(len(self.wf._data), 3)
        self.assertEqual(len(self.wf._data['x']), len(x))
        self.assertEqual(len(self.wf._data['y']), len(y))
        self.assertEqual(len(self.wf._data['e']), len(e))
        for j in range(len(x)):
            self.assertEqual(self.wf._data['x'][j], x[j])
            self.assertEqual(self.wf._data['y'][j], y[j])
            self.assertEqual(self.wf._data['e'][j], e[j])

    def test_execute(self):
        # indirectly tests normalise_grid and get_z_value

        # setup workflow + generate data
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_x_axis(0, 1, 2, 'x')
        self.wf.set_y_axis(1, 2, 2, 'y')
        self.func.add_function(ExpDecay())
        self.wf.set_scipy_engine([0, 0], [-9, -9], [9, 9])
        X, Y = self.wf.execute(self.func)

        grid = self.wf.get_grid
        expect_z = [[0.449, 1], [0, 0.115]]
        expect_x = [[0, 1], [0, 1]]
        expect_y = [[1, 1], [2, 2]]

        for i in range(2):
            for j in range(2):
                self.assertEqual(X[i][j], expect_x[i][j])
                self.assertEqual(Y[i][j], expect_y[i][j])
                self.assertAlmostEqual(grid[i][j],
                                       expect_z[i][j], 3)

    def test_get_slices(self):
        # setup workflow + generate data
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_x_axis(0, 1, 2, 'x')
        self.wf.set_y_axis(1, 2, 2, 'y')
        self.func.add_function(ExpDecay())
        self.wf.set_scipy_engine([0, 0], [-9, -9], [9, 9])
        _, _ = self.wf.execute(self.func)

        x, y = self.wf.get_slices()
        expect_x = [0.449, 1]
        expect_y = [1, 0.115]

        self.assertEqual(len(x), len(expect_x))
        self.assertEqual(len(y), len(expect_y))

        for j in range(len(x)):
            self.assertAlmostEqual(x[j], expect_x[j], 3)
            self.assertAlmostEqual(y[j], expect_y[j], 3)

    def test_get_parameters_and_errors(self):
        # rm
        pass

    def test_fails_if_no_data(self):
        # same
        with self.assertRaises(ValueError):
            self.wf.set_scipy_engine([0], [-9], [9])

    def test_execute_no_engine(self):
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_x_axis(0, 1, 2, 'x')
        self.wf.set_y_axis(1, 2, 2, 'y')
        with self.assertRaises(ValueError):
            _, _ = self.wf.execute(self.func)

    def test_execute_no_x_axis(self):
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_y_axis(1, 2, 2, 'y')
        with self.assertRaises(ValueError):
            _, _ = self.wf.execute(self.func)

    def test_execute_no_y_axis(self):
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_x_axis(0, 1, 2, 'x')
        with self.assertRaises(ValueError):
            _, _ = self.wf.execute(self.func)

    def test_add_second_engine_errors(self):
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_x_axis(0, 1, 2, 'x')
        self.wf.set_y_axis(1, 2, 2, 'y')
        self.wf.set_scipy_engine([], [], [])
        with self.assertRaises(RuntimeError):
            self.wf.set_scipy_engine([], [], [])

    def test_set_scipy_engine(self):
        # same
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([], [], [])
        self.assertEqual(self.wf.fit_engine.name, 'scipy')

    def test_update_scipy_fit_engine(self):
        # same
        self.assertEqual(self.wf.fit_engine, None)
        x, y, e = gen_grid_search_data()
        self.wf.preprocess_data(x, y, e)
        self.wf.set_scipy_engine([1.], [-4], [4])
        self.assertEqual(self.wf.fit_engine._guess, [1.])
        self.assertEqual(self.wf.fit_engine._lower, [-4])
        self.assertEqual(self.wf.fit_engine._upper, [4])

        bg = FlatBG()
        self.wf.update_fit_engine(bg, [2])
        self.assertEqual(self.wf.fit_engine._guess, [2.])
        self.assertEqual(self.wf.fit_engine._lower, [-4])
        self.assertEqual(self.wf.fit_engine._upper, [4])
