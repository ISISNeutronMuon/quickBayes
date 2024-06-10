Bayesian Methods
================

There are wide variety of Bayesian methods.
This section discusses the ideas and concepts of these other methods.
In addition it will highlight the key differences.



Bayes Theorm
------------

Bayesian inference is used to calculate the whole posterior probability distribution function (PDF).
The shape and spread of the PDF provides insight into the model and the parameter.
A narrow PDF suggests that the parameter is well defined or has little evidence of variation.
A broad PDF could mean that the model is insensitive to that specific parameter.
If the PDF is non-symmetric, then the most likely value is still the peak of the distribution.
However, it is more likely for the parameter to have a value higher/lower than the peak of the distribution.

Calculating the full PDF can be achieved by using Macov Chain Monte Carlo (MCMC) or nested sampling.
Essentially these methods will sample the PDF directly, allowing them to generate the full PDF.

Bayesian model selection use Bayes theorm to calculate the probability, :math:`P` of the data :math:`D` given the model :math:`M`

.. math::
    :name: eq:int

    P(D|M) = \int_\Omega P(D| \underline{\theta}, M)P( \underline{\theta}|M)\mathrm{d\underline{\theta}}.

where the :math:`\underline{\theta}` are the parameters and the integral is over all possible values for the parameters, :math:`\Omega`.

quickBayes
----------

The quickBayes method makes a series of assumptions to reduce :ref:`the full PDF evaluation <eq:int>` to a single analytic equation.
The full theory is discussed here.
The key assumptions are:

- The model can be written as a series of indistinguishable lines (i.e. the same repeated function)
- The lines can be written as an amplitude multiplied by some function
- The prior probabilities are flat across the domain of interest
- The normalisation is just the value of the hyper volume for the parameters.
- The :math:`P(D|\underline\theta, M)` can be represented as gaussians and be represented by a first order Taylor expansion


Odds factor
-----------

Check the label for :ref:`integration <eq:int>`

`odds factor<https://jakevdp.github.io/blog/2015/08/07/frequentism-and-bayesianism-5-model-selection/>`


MCMC to eval integral
---------------------

`MCMC <https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo>`

Nested sampling
---------------

`nested sampling <https://en.wikipedia.org/wiki/Nested_sampling_algorithm>`

AIC and BIC
-----------


`AIC <https://en.wikipedia.org/wiki/Akaike_information_criterion>`

`BIC <https://en.wikipedia.org/wiki/Bayesian_information_criterion>`
