Loglikelihood
=============

As discussed earlier the loglikelihood represents the probability of the model given the data.
To make the calculation faster a series of assumptions are made:

- The fit for the model is at a minima.
- The model can represent the data.

"Data Analysis A bayesian tutorial second edition", by D. S. Sivia (equation 4.20 page 88) uses the above assumptions to show that the probability can be written as

.. math::
    P(M|D_k,I) \propto \frac{N! (4\pi)^N}{\beta^N\sqrt{\det(H)}}\exp{(-\frac{\chi^2}{2})},

where :math:`\beta = (x_{max}-x_{min})A_{max}`, :math:`H` is the Hessian matrix.
The :math:`\beta` value represents the variation in the data by taking the product of the length of the :math:`x` range with the maximum height of the data.
To convert the probability into a loglikelihood we need to take the :math:`\log_{10}` and use:

.. math::
    \log_{10}(\exp{\alpha}) = \ln(\exp{\alpha})\log_{10}(\exp{1}) \\
    \log_{10}(\exp({\alpha}) = \alpha\log_{10}(\exp{1}) \\
    N! = \prod_{i=0}^{N} i => \log{N!} = \sum_{i=0}^N \log(i) \\
    \log(ab) = \log(a) + \log(b) \\
    \log(\frac{a}{b}) = \log(a) - \log(b) \\
    \log(a^N) = N\log(a)

to get:

.. math::
    \sum_{j=1}{N}\log(j) + N\log(4\pi) - 0.5\chi^2\log(\exp(1))
    - N\log(\beta) - 0.5\log(\det(H)).

This final equation is the unnormalized loglikelihood.

The assumption of being at a local minima is important, because if a fit is over-parameterised it should be unlikely (i.e. a large negative loglikelihood).
However, this is not the case as it produces a low :math:`\chi^2` value.
To help understand this imagine some data that averages around the value two.
If we fit a flat background to the data we would get the expected result of two.
If we then fit a sum of two flat background to the same data then the sum of the two will give the expected result of two.
However, this has added an extra degree of freedom as the values for the two flat backgrounds could be anything as long as they sum to two.
In terms of parameter space this represents a valley of possible values, hence it violates the assumption of being at a local minima.
To prevent this scenario a check is performed to check that the parameters are not too correlated and if they are a penalty is added to the loglikelihood.

Using the example from the introduction section we can calculate the loglikelihood for one and two gaussians.

.. code-block:: python

   import numpy as np
   from quickBayes.fit_functions.gaussian import Gaussian as func
   from quickBayes.functions.composite import CompositeFunction
   from quickBayes.fitting.scipy_engine import ScipyFitEngine
   from quickBayes.log_likelihood import loglikelihood
   # generate some data
   x = np.linspace(0,10, 100)
   noise = 1 + 0.1 * (np.random.normal(0, .2, len(x)))
   gauss = func()
   ground_truth = gauss (x, 103, 4.2, .9)
   y = ground_truth * noise
   e = np.sqrt(y) # normal for count data
   # construct fit functions
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


If the :code:`e` values are changed (e.g. :code:`e = np.power(y, 0.1)`) then the values for loglikelihoods will change.
It is even possible to change if one or two peaks are more likely given your data, by just altering the size of the error bars.
