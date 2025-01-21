.. _maths:

Bayesian theory
===============

In this section the concepts for using Bayesian techniques for model selection will be discussed.
The first part discusses Bayes theorm and how it can be used to get the probability of the model given the data.
The second part is on the odds factor and how to interpret the posterior probability for model selection.

Bayes Theorem
-------------

Bayesian inference is used to calculate the whole posterior probability distribution function (PDF).
The equation for the posterior probability can be written as

.. math::
   :name: bayes

   P(\underline{\theta} | D, M) = P(D | \underline{\theta}, M)\frac{P(D | \underline{\theta}, M)}{P(D | M)},

where :math:`\underline{\theta}` is a vector of model parameters, :math:`M` is the model and :math:`D` is the data.
:math:`P(\underline{\theta} | M)` is the prior distribution and represents current knowledge of the system.
:math:`P(D | M)` is the evidence and acts as a normalisation for the posterior.
The shape and spread of the PDF provides insight into the model and the parameter.
A narrow PDF suggests that the parameter is well defined or has little evidence of variation.
A broad PDF could mean that the model is insensitive to that specific parameter.
If the PDF is non-symmetric, then the most likely value is still the peak of the distribution.
However, it is more likely for the parameter to have a value higher/lower than the peak of the distribution.

Calculating the full PDF can be achieved by using :ref:`Macov Chain Monte Carlo (MCMC) <MCMC>` or :ref:`nested sampling<nest>`.
Essentially these methods will sample the PDF directly, allowing them to generate the full PDF.

Bayesian model selection use Bayes theorm to calculate the probability, :math:`P` of the data :math:`D` given the model :math:`M`

.. math::
   :name: eq_int

   P(D|M) = \int_\Omega P(D| \underline{\theta}, M)P( \underline{\theta}|M)\mathrm{d\underline{\theta}}.

where the :math:`\underline{\theta}` are the parameters and the integral is over all possible values for the parameters, :math:`\Omega`.

Odds factor
-----------

One method for comparing models is known as the odds factor.
It assumes that you can calculate the probability of the data given the :math:`i^{\mathrm{th}}` model (:math:`M_i`), by taking the ratio of two different models.
This section will use a derivation based on the work from `here <https://jakevdp.github.io/blog/2015/08/07/frequentism-and-bayesianism-5-model-selection/>`_.

For model selection we want the model posterior

.. math::
   P(M | D) = P(D | M) \frac{P(M)}{P(D)},

where :math:`P(D | M)` is the probability of the data given the model, :math:`P(M)` is the probability of the model and :math:`P(D)` the probability of the data.
The probability of the data will be the same for all models, so by taking a ratio it cancels out

.. math::
   :label: odds

   O_{21} = \frac{P(M_2 | D)}{P(M_1 | D)} = \frac{P(D | M_2)P(M_2)}{P(D | M_1)P(M_1)}

where :math:`0_{21}` is the odds factor for models two (:math:`M_2`) and one (:math:`M_1`).
Assuming that there is no prior knowledge then :math:`P(M_1) \approx P(M_2)`.
Then equation :math:numref:`odds` can be simplified to

.. math::
   O_{21} = \frac{P(D | M_2)}{P(D | M_1)},

which is known as the Bayes factor.
Alternatively, the Bayesian probability for the :math:`j^\mathrm{th}` model can be written as

.. math::
   P(M_j | D) = \frac{ P(D | M_j)}{ \sum_k P(D | M_k)}.


To evaluate the odds factor, the probability of the data given the model needs to be calculated.
This is typically written as

.. math::
   :label: P_int

   P(D | M) = \int_\Omega d\underline{\theta} \quad P(D| \underline{\theta}, M)P(\underline{\theta} | M)

where the integral over :math:`\Omega` is over the available parameter space for :math:`\underline{\theta}`.
This quantity can be evaluated directly using either :ref:`Markov Chain Monte Carlo (MCMC) <MCMC>` or :ref:`nested sampling <nest>`.

