Fit Functions
=============

Fit functions are the mathematical representations of the data.
In quickBayes the fit functions are all classes with a very specific structure, see :py:mod:`here<quickBayes.fit_functions>` for a full list of currently available options.
The first step is to create a fit function object, these are Python classes designed to make it easier for inspecting important information.
As an example we will use the gaussian function and this page will outline the various methods that you can use.
The simplest example is to just create a fit function and to then call it.

.. code-block:: python

  from quickBayes.fit_functions.gaussian import Gaussian
  import numpy as np

  function = Gaussian()
  x = np.linspace(-5, 5, 100)
  # call the function
  y = function(x, amplitude=1, x0=0, sigma=0.2)

This provides an easy way to evaluate the function.
You can also check the number of parameters by using :code:`function.N_params` to check.
However, we may want to keep the current parameters for later (e.g. they are the results of a fit).
This can be done by;

.. code-block:: python

  results = {}
  results = function.report(results, 1, 0, .2)
  print(results)

which will output

.. code-block:: python

  {'Amplitude': [1.0], 'Mean': [0], 'Sigma': [2.0]}

Repeated calls to the report method will append to the lists

.. code-block:: python

  results = function.report(results, 2, 1, .2)
  print(results)

Will output:

.. code-block:: python

  {'Amplitude': [1.0, 2], 'Mean': [0, 1], 'Sigma': [2.0, 0.2]}

It is also possible to then read the results from the report output.
This makes it easy to extract the parameters from the report (e.g. after multiple fits).
By default the method will read the first entry into the report

.. code-block:: python

  params = function.read_from_report(report)
  print(params)

Hence,

.. code-block:: python

  [1.0, 0.0, 2.0]

Or a specific index can be given

.. code-block:: python

  params = function.read_from_report(report, 1)
  print(params)

which gives

.. code-block:: python

  [2.0, 1.0, 0.2]

Similar to the `report` method is a `report_errors` method:

.. code-block:: python

   errors = [ .1, .01, .02]
   error_report = function.report_errors({}, errors, params)

The behaviour for this function is the same as :code:`report`, but some functions have non-trivial errors that will be calculated as part of this method.

When doing a fit, it can be useful to get an initial guess value :code:`guess = function.get_guess()`.
If the guess is not appropritate it can be changed by using :code:`function.set_guess([10, 4, 1])`.
Similarly the bounds for the function are given by :code:`lower, upper = function.get_bounds()`.
The bounds can be set by using :code:`function.set_bounds([1, 2, 0.4], [10, 7, 1])`.

