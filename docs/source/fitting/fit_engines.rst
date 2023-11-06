Fit Engines
===========

The fit engine provides a wrapper to allow for a harmonized interface to different optimization methods.
Making it easier to swap the optimization method in a script to a different one.
To make it easier to compare to real data the fit engine is initialized with the :math:`x, y, e` data that you are interested in (e.g. experimental data).
When doing a fit the :math:`x’, y’, e’` data are provided and can be different (e.g. rebinned) and these are the values that the fit is performed against.
However, when the results are reported a spline is used to map the fit back onto the original :math:`x` axis.
A single fit engine instance can be used to calculate multiple fits and it will remember the full history.
All of the access methods have an :code:`index` argument that allows access to the history (note that it starts with :math:`0`).
The access methods are:

- :code:`Get_chi_squared`.
- :code:`Get_covariance_matrix`.
- :code:`Get_fit_values` returns the :math:`x` data, the splined fit values, the splined fit error bars, the difference between :math:`y` and the fit and the error of the difference

The fit is called using :code:`do_fit` which takes :math:`x, y, e` data and a fit function object.

At present there are two fit engines, :code:`ScipyFitEngine` and :code:`GoFitEngine`.


Scipy
=====

The :code:`ScipyFitEngine` uses the :code:`curve_fit` method from scipy.
At initialisation it also requires

- The lower bounds for the fit parameters.
- The upper bounds for the fit parameters.
- The initial guess for the fit parameters.
- The maximum number of iterations (defaults to :math:`220000`).

There is also a :code:`set_guess_and_bounds` method to update the guess and bounds values.

.. code-block:: python

    import numpy as np
    from quickBayes.fit_functions.gaussian import Gaussian as func
    from quickBayes.fitting.scipy_engine import ScipyFitEngine
    	# generate some fake data
    x = np.linspace(0,10, 100)
    noise = 1 + 0.1*(np.random.normal(0, .2, len(x)))
    g_func = func()
    ground_truth = g_func(x, 103, 4.2, .9)
    y = ground_truth*noise
    e = np.sqrt(y) # errors for count data

    g_func2 = func()
    # do the fit
    engine = ScipyFitEngine(x, y, e, lower=[0, 0, 0], upper=[200, 10, 10], guess=[100, 5, 1.2])
    engine.do_fit(x, y, e, g_func2)
    chi_2 = engine.get_chi_squared()
    params = engine.get_fit_parameters()
    print(chi_2, params)


GoFit
=====

The :code:`GoFitEngine` uses the :code:`multistart` method from `GoFit <https://ralna.github.io/GOFit/_build/html/index.html>`.
At initialisation it also requires

- The lower bounds for the fit parameters.
- The upper bounds for the fit parameters.
- The number of samples (defaults to :math:`10`).
- The maximum number of iterations (defaults to :math:`220000`).

It also has a :code:`set_bounds_and_N_params` method for updating the bounds and the number of samples.

