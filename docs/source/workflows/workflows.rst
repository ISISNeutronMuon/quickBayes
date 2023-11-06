Workflow
========

In the loglikelihood section we calculated if one or two gaussians were more likely for a dataset.
To make user scripts easier to work with quickBayes has implemented a :code:`workflow` class.
This class allows for the complex operations required by an analysis to be implemented once.
When its imported into a user script the workflow significantly simplifies the code making it easier to read and understand.
For example lets consider the case of wanting to know how many gaussians are within a dataset.
We know that there is at least one gaussian and no more than five.
All of the peaks are centred near zero and are approximately zero outside of the range :math:`-1` to :math:`1`.
The x data starts at :math:`-5` and ends at :math:`5`, but most of this is just background.
Hence, it is worth cropping the data before the main part of the analysis.
Assuming the workflow exists (given in the next section) the following could be used:


.. code-block:: python

   import numpy as np
   from quickBayes.fit_functions.gaussian import Gaussian as func
   from quickBayes.functions.composite import CompositeFunction
   from my_workflow import GaussianExample

   # generate some data
   x = np.linspace(-50, 50, 1000)
   noise = 1 + 0.1 * (np.random.normal(0, .2, len(x)))
   gauss = Gaussian()
   ground_truth = gauss (x, 103, 4.2, .9)
   y = ground_truth * noise
   e = np.power(y, 0.5) # normal for count data
   
   # setup 
   params = {}
   errors = {}
   func = CompositeFunction() # could add a background in needed
   workflow = GaussianExample(params, errors)
   # this will crop the data
   workflow.preprocess_data(x, y, e, -10., 10.)
   workflow.set_scipy_engine([], [], [])
   # execute the workflow with no intial guess and up to 5 gaussians
   _ = workflow.execute(max_num_features=5, func=func, params=[])
   
   # report results
   params, errors = workflow.get_parameters_and_errors
   print(params.keys())
   for j in range(5):
   
       print(f"{j + 1} peak(s)")
       print('loglikelihood', params[f"N{j+1}:loglikelihood"])
       for key in errors.keys():
           if f"N{j+1}" in key:
               print(key, params[key], errors[key])

Every workflow inherits the :code:`workflowTemplate` class.
This has a :code:`fit_engine` property, which is set with one of the following commands:

- :code:`set_scipy_engine`.
- :code:`set_gofit_engine`.

The other methods are

- :code:`preprocess_data` for any preprocessing that might be required (e.g. cropping or rebinning the data).
- :code:`update_fit_engine` allows for updates to be passed to the specified fit engine (e.g. a new guess for the scipy fit engine).
- :code:`execute` for doing the analysis.


Model Selection
===============

The main advantage of a workflow is to simplify model selection.
Examples for the models being selected between could be the number of Lorentzians (:code:`QLData`) or the number of exponential decays (:code:`MuonExpDecay`) .
The :code:`ModelSelectionWorkflow` simplifies creating these kind of analyses.
It inherits the the :code:`workflowTemplate` class but has as additional property :code:`get_parameters_and_errors`.
It also has two additional methods

- :code:`update_function` for updating the model (e.g. adding a peak).
- :code:`report` for updating dictionaries with the results.

For example lets consider the case of wanting to know how many gaussians are within a dataset.
We know that there is at least one gaussian.
All of the peaks are centred near zero and are approximately zero outside of the range :math:`-10` to :math:`10`.
The x data starts at :math:`-50` and ends at :math:`50`, but most of this is just background.
Hence, it is worth cropping the data before the main part of the analysis.
The following workflow could be used:

.. code-block:: python

    import numpy as np
    from quickBayes.fit_functions.gaussian import Gaussian
    from quickBayes.utils.crop_data import crop
    from quickBayes.workflow.model_template import ModelSelectionWorkflow
    from quickBayes.functions.base import BaseFitFunction
    from numpy import ndarray
    

    class GaussianExample(ModelSelectionWorkflow):
        """
        A class for the finding gaussians
        """
        def preprocess_data(self, x_data: ndarray,
                            y_data: ndarray, e_data: ndarray,
                            start_x: float, end_x: float) -> None:
            """
            The preprocessing needed for the data.
            This crops and stores the data.
            :param x_data: the x data to fit to
            :param y_data: the y data to fit to
            :param e_data: the errors for the y data
            :param start_x: the start x value
            :param end_x: the end x value
            """
            sx, sy, se = crop(x_data, y_data, e_data,
                              start_x, end_x)
            super().preprocess_data(sx, sy, se)
    
        @staticmethod
        def _update_function(func: BaseFitFunction) -> BaseFitFunction:
            """
            This method adds a Gaussian to the fitting
            function.
            :param func: the fitting function that needs modifying
            :return the modified fitting function
            """
    
            g_function = Gaussian()
            # need to change the bounds and guess on the function
            g_function.set_bounds([80, 0, .2], [120, 10, 2])
            g_function.set_guess([100, 5, 1.2])
            func.add_function(g_function)
            return func

The following workflows are available as part of the quickBayes package:

- :code:`QLData` for determining if 1, 2 or 3 Lorentzians are present in qausielastic data.
- :code:`QlStretchedExp` for getting for loglikelihood of a single stretched exponential for quasielastic data.
- :code:`MuonExpDecay` for determining if 1, 2, 3 or 4 decays are present in MuSR data.

All of these workflows use the scipy fit engine.


Grid Search
===========

It is possible to use a workflow to do a grid search of two parameters.
This requires a function with at least one free parameter, after fixing two for the grid point. 
At present this is an experimental method and is not fully supported. 
The fitting function needs to be made specifically for this method, with the fixed paramters as private members. 
The results from the grid search are the loglikelihood normalised so that the most likely point has a value of one.
The least likely point has a value of zero.
It is important to note that this normalisation only takes into account the sampled grid.
Hence, if the grid is too coarse then the true most likely value will not be found. 
