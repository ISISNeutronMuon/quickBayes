.. _MCMC:

Markov Chain Monte Carlo (MCMC)
-------------------------------

The integral for the evidence is

.. math::

   P(D|M) = \int_\Omega P(D| \underline{\theta}, M)P( \underline{\theta}|M)\mathrm{d\underline{\theta}}

which typically requires a numerical method for its evaluation.

Markov Chain Monte Carlo (MCMC) is a method that uses random walkers to estimate the probability distribution.
The MCMC has two main components, the first defines how the walkers select their new positions and the second is how to determine if to accept the new values.
There are several options for each of these parts, leading to numerous possible MCMC simulations.
For this section, the differential evolution and Metropolis Hastings methods will be considered.
These methods are relatively straightforward and provide insight into MCMC methods.
More detailed discussion can be found `here <https://www.sciencedirect.com/science/article/pii/S0022169423007643>`_ and `here <https://royalsocietypublishing.org/doi/full/10.1098/rsta.2014.0405>`_.

The differential evolution algorithm provides an equation that determines how the random walkers will be updated.
Each walker describes the potential parameters for the model, :math:`\underline{\theta}_i^t`, where :math:`i` is used to label the different walkers and :math:`t` labels the step of the algorithm.
The walkers evolve to their proposed new positions according to the equation

.. math::
   \underline{\theta}_i^{t+1} = \underline{\theta}_i^t + \gamma(\underline{\theta}_j^t - \underline{\theta}_k^t) + \underline{\epsilon}

where :math:`i\ne j\ne k`, :math:`\gamma` gives the strength of the coupling between the walkers and :math:`\epsilon` provides a random change to the parameters.
The first term provides some history to the method, so the new parameter values are not completely random.
This prevents the walker from picking new guesses that are significantly worse than the current one.
The second term can be thought of as a diffusion term, it determines how much the walker will move.
It depends on the distance between two different walkers and then multiplies the result by a scalar :math:`\gamma`.
So if the MCMC has a bad estimate of the PDF, then the distance between the walkers is probably large.
Hence, the next guess will move more in an attempt to find a better set of parameters.
This part of the algorithm is described as the burn in period, and is the time for the walkers to find a good PDF.
If the walkers provide a good description of the PDF, then the difference is small.
This results in the walker moving a small distance from its original position.
The competing effects from the second term, depending on the distance between two walkers, pulls the walkers towards good parameter values.
The final term adds a bit of randomness into the position update.

Once a walker has a new position, it will not automatically move to it.
Instead it has a finite probability of its values being updated.
One method for determining if the walker position should be updated is the Metropolis Hastings algorithm.
For a given walker position it is possible to calculate the likelihood (typically a Gausssian).
The likelihoods are calculated for both the proposed, :math:`\underline{\theta}_j^{t+1}`, and current, :math:`\underline{\theta}_j^t`, walker positions.
Then the ratio is taken and compared to a random number.
The proposed value is accepted if the random number is smaller than the ratio.
As a result the walker is not guranteed to update to the new position if its current values are good enough.

The combined effect of the two algorithms is that the walkers will eventually converge onto a distribution that describes the PDF.
To make this tractable, a few extra conditions are needed.
The first is to limit the probability space, by placing limits on the potential parameter values.
This is the prior for the problem and the starting distribution for the walkers is normally flat across the parameter space.
The second is to define the behaviour of the walkers at the boundaries of the parameter space.
Typically they are chosen to be either reflective or a periodic boundary.
For a Mertropolis Hastings algorithm both options are suitable, because the results are independent of the path taken by the walker.

MCMC will eventually give a good representation of a unimodal posterior PDF of the data, even if it has some complex structure.
However, it can be very computationally expensive to evaluate due to the number of walkers required to get a good estimate of the posterior and the burn in period.
The compuational cost is difficult to estimate as it requires prior knowledge of how long the burn in period should be.
Hence, it can be too short leading to poor results or too long wasting valuable computational time.


