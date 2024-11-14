.. _cf:

Bayesian Methods
================

There are wide variety of Bayesian methods.
This section discusses the ideas and concepts of these other methods.
In addition it will highlight the key differences.



Bayes Theorm
------------

Bayesian inference is used to calculate the whole posterior probability distribution function (PDF).
The equation for the posterior probability can be written as

.. math::
   P(\underline{\theta} | D, M) = P(D | \underline{\theta}, M)\frac{P(D | \underline{\theta}, M)}{P(D | M)},

where :math:`\underline{\theta}` is a vector of model parameters, :math:`M` is the model and :math:`D` is the data.
:math:`P(\underline{\theta} | M)` is the prior distribution and represents current knowledge of the system.
:math:`P(D | M)` is the evidence and acts as a normalisation for the posterior.
The shape and spread of the PDF provides insight into the model and the parameter.
A narrow PDF suggests that the parameter is well defined or has little evidence of variation.
A broad PDF could mean that the model is insensitive to that specific parameter.
If the PDF is non-symmetric, then the most likely value is still the peak of the distribution.
However, it is more likely for the parameter to have a value higher/lower than the peak of the distribution.

Calculating the full PDF can be achieved by using Macov Chain Monte Carlo (MCMC) or nested sampling.
Essentially these methods will sample the PDF directly, allowing them to generate the full PDF.

Bayesian model selection use Bayes theorm to calculate the probability, :math:`P` of the data :math:`D` given the model :math:`M`

.. math::
   :name: eq_int

   P(D|M) = \int_\Omega P(D| \underline{\theta}, M)P( \underline{\theta}|M)\mathrm{d\underline{\theta}}.

where the :math:`\underline{\theta}` are the parameters and the integral is over all possible values for the parameters, :math:`\Omega`.

quickBayes
----------

The quickBayes method makes a series of assumptions to reduce :ref:`the full PDF evaluation <eq_int>` to a single analytic equation.
The full theory is discussed here.
The key assumptions are:

- The model can be written as a series of indistinguishable lines (i.e. the same repeated function)
- The lines can be written as an amplitude multiplied by some function
- The prior probabilities are flat across the domain of interest
- The normalisation is just the value of the hyper volume for the parameters.
- The :math:`P(D|\underline\theta, M)` can be represented as gaussians and be represented by a first order Taylor expansion


Odds factor
-----------

One method for comparing models is known as the odds factor.
It assumes that you can calculate the probability of the data given the :math:`i^{\mathrm{th}` model (:math:`M_i`), by taking the ratio of two different models.
This section will use a derivation based on the work from `here <https://jakevdp.github.io/blog/2015/08/07/frequentism-and-bayesianism-5-model-selection/>`_.

For model selection we want the model posterior

.. math::
   P(M | D) = P(D | M) \frac{P(M)}{P(D)},

where :math:`P(D | M)` is the probability of the data given the model, :math:`P(M)` is the probability of the model and :math:`P(D)` the probability of the data.
The probability of the data will be the same for all models, so by taking a ratio the term can be removed

.. math::
   :label: odds

   O_{21} = \frac{P(M_2 | D)}{P(M_1 | D)} = \frac{P(D | M_2)P(M_2)}{P(D | M_1)P(M_1)}

where :math:`0_{21}` is the odds factor for models two (:math:`M_2`) and one (:math:`M_1`).
Assuming that there is no prior knowledge then :math:`P(M_1) \approx P(M_2)`.
Then equation :math:numref:`odds` can be simplified to

.. math::
   O_{21} = \frac{P(D | M_2)}{P(D | M_1)},

which is known as the Bayes factor.
Alternatively, the Bayesian probability for the :math:`j^\mathrm{th}` model is

.. math::
   P(M_j | D) = \frac{ P((D | M_j)}{ \sum_k P(D | M_k)}.


To evaluate the odds factor, the probability of the data given the model needs to be calculated.
This is written as

.. math::
   :label: P_int
   P(D | M) = \int_\Omega d\underline{\theta} \quad P(D| \underline{\theta}, M)P(\underline{\theta} | M)

where the integral over :math:`\Omega` is over the available parameter space for :math:`\underline{\theta}`.
This quantity can be evaluated using either Markov Chain Monte Carlo (MCMC) or nested sampling.



Makov Chain Monte Carlo (MCMC)
------------------------------

The integral in equation :math:numref:'P_int' typically requires a numberical method to evaluate it.
Markov Chain Monte Carlo (MCMC) is a method that uses random walkers to estimate the probability distribution.
The MCMC has two main components, the first defines how the walkers select their new positions and the second is how to determine if to accept the new values.
There are several options for each of these parts, leading to numerous possible MCMC simulations.
For this section, the differential evolution and Metropolis Hastings methods will be considered.
These methods are relatively straightforward and provide insight into MCMC methods.

The differential evolution algorithm provides an equation that determines how the random walkers will be updated.
Each walker describes the potential parameters for the model, :math:`\underline{\theta}_i^t`, where :math:`i` is used to label the different walkers and :math:`t` labels the step of the algorithm.
The walkers are then evolve to their new positions according to the equation

.. math::
   :label:MCMCDE
   \underline{\theta}_i^{t+1} = \underline{\theta}_i^t + \gamma(\underline{\theta}_j^t - \underline{\theta}_k^t) + \underline{\epsilon}

where :math:`i\ne j\ne k`, `gamma` gives the strength of the coupling between the walkers and :math:`\epsilon` provides a random change to the parameters.
The first term provides some history to the method, so the new parameter values are not completely random.
This prevents the walker from picking new guesses that are significantly worse than the current one.
The second term can be thought of as a diffusion term, it determines how much the walker will move.
It depends on the distance between two different walkers and then multiplies the result by a scalar :math:`gamma`.
So if the MCMC has a bad estimate of the PDF, then the distance between the walkers is probably large.
Hence, the next guess will move more in an attempt to find a better set of parameters.
This part of the algorithm is described as the burn in period, and is the time for the walkers to find a good PDF.
If the walkers provide a good description of the PDF, then the difference is small.
This results in the walker moving very far from its original position.
The competing effects from the second term, depending on the distance between two walkers, has the effect of pulling walkers towards good parameter values.
The final term just randoms a bit of randomness into the position update.

Once a walker has a new position, it will not automatically move to it.
Instead it has a finite probability of its values being updated.
One method for determining if the walker position should be updated is the Metropolis Hastings algorithm.
For a given walker position it is possible to calculate the likelihood (typically a Gausssian).
The likelihoods are calculated for both the proposed new (:math:`\underline{\theta}_j^{t+1}`) and current (:math:`\underline{\theta}_j^t`) walker positions.
Then the ratio is taken and compared to a random number.
The proposed value is accepted if the random number is smaller than the ratio.
As a result the walker is not guranteed to update to the new position if its current values are good.

The combined effect of the two algorithms is that the walkers will eventually converge onto a distribution that describes the PDF.
To make this tractable, a few extra conditions are needed.
The first is to limit the probability space, by placing limits on the potential parameter values.
This is the prior for the problem and the starting distribution for the walkers is normally flat across the parameter space.
The second is to define the behaviour of the walkers at the boundaries of the parameter space.
Typically they are chosen to be either reflective or a periodic boundary.
For a Mertropolis Hastings algorithm both options are suitable, because the results are independent of the path taken by the walker.

MCMC will eventually give a very good represnetation of the data, even if it has a complex PDF.
However, it can be very computationally expensive to evaluate due to the burn in period.
It is also difficult to know prior to the calculation how long the burn in period should be.
Hence, it can be too short leading to poor results or too long wasting valuable computational time.


`MCMC <https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo/>`_

Nested sampling
---------------

`nested sampling <https://en.wikipedia.org/wiki/Nested_sampling_algorithm/>`_

AIC and BIC
-----------

`AIC <https://en.wikipedia.org/wiki/Akaike_information_criterion/>`_

`BIC <https://en.wikipedia.org/wiki/Bayesian_information_criterion/>`_
