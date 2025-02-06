.. _stats:

Introduction
============

Statistics is dominated by a school of thought known as frequentist statistics.
To demonstrate the difference, let's consider a bag containing a mixture of red and yellow balls.
If a ball is chosen from the bag, such that only one person knows its colour.
Then a frequentist would say that the probability of that ball being red is either :math:`100\%` or :math:`0\%`, as it is definitely red or not.

A Bayesian approach relies more on probability and uncertainty.
In the above example of picking a ball, they would say that the probability of the ball being red is :math:`50\%` as it is either red or yellow.
This initial probability is known as a prior and represents the known information at the beginning.
After revealing the colour of the ball it is placed back into the bag, it is shaken, and a new ball is picked.
After repeating this :math:`100` times a red ball was picked :math:`70` times.
This leads to a hypothesis that :math:`70%` of the balls are red.
After another :math:`100` draws this hypothesis can be tested by using Bayes theorem to calculate the probability of the hypothesis given the new data.
This probability is known as the posterior probability and represents the likelihood of a hypothesis given the data.
Bayes theorem for the posterior probability is

.. math::

   P(\theta|y) = \frac{P(y|\theta)P(\theta)}{P(y)},

where :math:`\theta` are the model paramters and :math:`y` are the measured data.


Why use Bayes?
--------------

In science it is often difficult to know which model gives the best representation of the measured data.
For example, in muon spectroscopy it is important to know how many exponential decays are present in the data, with each one providing information on the muon stopping site or the muon relaxation time (depending on the setup).
The model is used to describe the hypothesis, in this case there are :math:`M` exponential decays, where :math:`M` is an integer between :math:`1` and :math:`N`.
This results in :math:`N` different models and knowing which one to use is not obvious.
Some methods for deciding include:

- Predictions from theory
- Best chi squared value
- Bayesian methods

However, even theory will sometimes be unable to identify a single correct model.
Instead, it may limit the options, for example its one or two decays, and something else needs to be used to determine the correct model.
Using non-linear least squares to calculate the chi squared value for the model and data, then naively selecting the model with the lowest value is another option.
This can be misleading if not done carefully as the model may be over-parameterized, resulting in physically meaningless parameters but a good chi squared value.
The last option is to use a Bayesian method, which calculates the probability of the model given the data.
It does not state that a given model is correct, just that it is more likely.
