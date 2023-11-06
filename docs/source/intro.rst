Introduction
============

The quickBayes package is an open source library for calculating Bayesian quantities in a short period of time, by making some assumptions. 
The package is cross platform, supporting Windows, Mac OS and Linux. 
This package has been developed by Anthony Lim from STFC’s ISIS Neutron and Muon facility. 

Background
==========

A Bayesian analysis answers the question; what is the probability of my model given my data?
The first part is easy to understand, if the model is a linear function and the data is a peak then the model is unlikely to be a good choice. 
The second part is more subtle, the important part is; given my data.
To help explain this we will consider two sets of data, they will have identical x and y values and the only difference will be the errors associated with the y values. 
For this example we will use a single Gaussian with an amplitude of :math:`103`, mean of :math:`4.2` and sigma of :math:`0.9`.
The following code was used to produce the data

.. code-block:: python

  from quickBayes.fit_functions.gaussian import Gaussin
  import numpy as np

  x = np.linspace(0, 10, 100)
  noise = 1 + 0.1*(np.random,normal(0, 2, len(x)))
  function = Gaussian()
  ground_truth = function(x, 103, 4.2, .9)
  y = ground_truth * noise
  e1 = np.sqrt(y)
  e2 = np.power(y, 0.1)
     
The two data sets are shown below. 

.. figure:: /images/bayes_example_data.png
   :alt: bayes_example_data.png
   :width: 400px
   :align: center

Both data sets return the exact same fitting parameters with values :math:`102. \pm 3.` for the amplitude, :math:`4.20 \pm 0.03` for the mean and :math:`0.90 \pm 0.02` for sigma.
The chi squared values are 0.003 for the first data set (:math:`\sqrt(y)`) and 0.04 for the second data set (with the smaller error bars).
This is expected, but both results suggest a good fit. 

The loglikelihood
=================

Alternatively, we can use the loglikelihood to calculate the log of the probability that the fit describes the data.
Since it is a log of a probability the values will always be negative, with the most likely outcomes being closest to zero. 
For the first data set the loglikelihood is :math:`-275.8` and for the second data set it is :math:`-7.0`.
This means that the second data set is more likely to be represented by the model (fit), which is a consequence of the smaller error bars. 
Furthermore, if we then consider the possibility of two gaussians the loglikelihoods are :math:`-58.2` for the first data set and :math:`-59.4` for the second data set. 
These values are similar, but it is important to notice that for the first data set it is more likely to be two gaussians than one. 
Whereas the second data set is more likely to be one gaussian, as expected.
Traditionally Bayesian methods use some sort of Monte Carlo to sample the probability space for each parameter of interest from which it can calculate the loglikelihood. 
This is computationally very slow and can give the wrong results if the calculation has not had sufficient time to pass the burn-in period. 
The quickBayes package makes some assumptions to remove the need for a Monte Carlo calculation by making the equation for the loglikelihood analytic. 

In the following documentation we will cover the important parts of quickBayes using this example and how it can be used for model selection (e.g. for determining the number of peaks in some data).

