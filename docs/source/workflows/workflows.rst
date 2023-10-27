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
   from my_workflow import GaussianExample

   # generate some data
   x = np.linspace(-5, 5, 1000)
   noise = 1 + 0.1 * (np.random.normal(0, .2, len(x)))
   gauss = func()
   ground_truth = gauss (x, 103, 4.2, .9)
   y = ground_truth * noise
   e = np.sqrt(y) # normal for count data

   #
   workflow = GaussianExample(

   gauss_2 = func()
   two_peaks = CompositeFunction()
   two_peaks.add_function(gauss)
   two_peaks.add_function(gauss_2)
   # do fitting
   engine = ScipyFitEngine(x, y, e, lower=[0, 0, 0], upper=[200, 10, 10], guess=[100, 5, 1.2])
   engine.do_fit(x, y, e, g)
   chi_2 = engine.get_chi_squared()
   # calculate loglikelihood
   beta = np.max(y) * (np.max(x) - np.min(x))
   loglike = loglikelihood(len(x), chi_2, engine.get_covariance_matrix(), 1, beta)
   # report results
   print("one peak", engine.get_fit_parameters(), chi_2, loglike)

   # for two peaks
   engine.set_guess_and_bounds([100, 5, 1.2, 100, 5, 1.2], [0, 0, 0, 0, 0, 0], [200, 10, 10, 200, 10, 10])
   engine.do_fit(x, y, e, two_peaks)
   chi_2 = engine.get_chi_squared()
   loglike = loglikelihood(len(x), chi_2, engine.get_covariance_matrix(), 1, beta)
   print("two peaks", engine.get_fit_parameters(), chi_2, loglike)




Every workflow inherits the :code:`workflowTemplate` class.
This has a :code:`fit_engine` property, which is set with one of the following commands:

- :code:`set_scipy_engine`.
- :code:`set_gofit_engine`.

The other methods are

- :code:`preprocess_data` for any preprocessing that might be required (e.g. cropping or rebinning the data).
- :code:`update_fit_engine` allows for updates to be passed to the specified fit engine (e.g. a new guess for the scipy fit engine).
- :code:`execute` for doing the analysis.


Model Selection
=============

The main advantage of a workflow is to simplify model selection.
Examples for the models being selected between could be the number of Lorentzians (:code:`QLData`) or the number of exponential decays (:code:`MuonExpDecay`) .
The :code:`ModelSelectionWorkflow` simplifies creating these kind of analyses.
It inherits the the :code:`workflowTemplate` class but has as additional property :code`get_parameters_and_errors`.
It also has two additional methods

- :code:`update_function` for updating the model (e.g. adding a peak).
- :code:`report` for updating dictionaries with the results.

For example lets consider the case of wanting to know how many gaussians are within a dataset.
We know that there is at least one gaussian.
All of the peaks are centred near zero and are approximately zero outside of the range :math:`-10` to :amth:`10`.
The x data starts at :math:`-50` and ends at :math:`50`, but most of this is just background.
Hence, it is worth cropping the data before the main part of the analysis.
The following workflow could be used:

.. code-block:: python
    from quickBayes.functions.composite import CompositeFunction
    from quickBayes.functions.gaussian import Gaussian
    from quickBayes.utils.crop_data import crop
    from quickBayes.workflow.model_template import ModelSelectionWorkflow
    from quickBayes.functions.base import BaseFitFunction
    from numpy import ndarray
    from typing import Dict, List


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
            func.add_function(g_function)
            return func

The following workflows are available as part of the quickBayes package:

- :code:`QLData` for determining if 1, 2 or 3 Lorentzians are present in qausielastic data.
- :code:`QlStretchedExp` for getting for loglikelihood of a single stretched exponential for quasielastic data.
- :code:`MuonExpDecay` for determining if 1, 2, 3 or 4 decays are present in MuSR data.

All of these workflows use the scipy fit engine.

