.. _nest:

Nested sampling
---------------

A popular alternative to MCMC is nested sampling.
The algorithm creates a set of randomly distributed samples across the potential parameter space, just like MCMC.
However, the samples (refered to as walkers in MCMC) do not evolve in nested sampling.
Instead they are used to create a series of contours of approximatly equal likelihood within the parameter space.
This can be thought of as being similar to Russian dolls, where the larger outer shells are removed to reveal a smaller more complex shell.
As a result nested sampling is good for investigating multi-modal posterior distributions.

This is a brief description of how the algorithm works, but a more detailed discussion of the subject is outlined in this `paper <https://arxiv.org/pdf/2205.15570>`_.

The likelihood, :math:`P(D|\underline{\theta}, M)`, and prior, :math:`P(\underline{\theta}| M)`, are related to the evidence by

.. math::
   P(D|M) = \int_\Omega P(D| \underline{\theta}, M)P( \underline{\theta}|M)\mathrm{d\underline{\theta}}.

In nested sampling the notation is slightly different:

* The evidence is :math:`Z = P(D | M)`.
* The prior is :math:`\pi(\underline{\theta}) = P(\underline{\theta} | M)`.
* The liklihood is :math:`L(\underline{\theta}) = P(D | \underline{\theta}, M)`.

Hence, the above equation is written as

.. math::
   :label: NS

   Z = \int_\Omega L(\underline{\theta})\pi(\underline{\theta}) \mathrm{d\underline{\theta}},

which can be simplified to a one dimensional integral with a change of variables.
We define the likelihood contour to be

.. math::
   :label: contour

   X(L) = \int \pi(\underline{\theta}) \mathrm{d\underline{\theta}},

and the integral is over a surface with a constant likelihood.
In practice the integral needs to be replaced by a summation and then :math:`X` is described as the volume variable.
The change of variables allows equation :math:numref:`NS` to be written as

.. math::
   :label: Z

   Z = \int_0^1 L(X) \mathrm{dX}.

The first few steps are similar to MCMC:

#. The bounds for the parameter space are defined.
#. A set of uniformally random points are placed in the bound parameter space.
#. The likelihoods are calculated for all of the random points.

Nested sampling then simplifies the problem of exploring multidimensional space, by reducing it to a series of shells (contours in the limit of infinite samples, see equation :math:numref:`contour`).
This is done by initialising the evidence to zero and the volume variable to one.
The following set of steps are then repeated until a stopping criteria is met:

#. The sample with the minimum value for the likelihood, :math:`L*` is identified for an iteration :math:`i`.
#. The integral in equation :math:numref:`Z` is updated with the new likelihood, via a numberical integration method. For trapezium rule the new contribution to the evidence will be :math:`\frac{L*(X_{i-1} - X_{i+1})}{2}`.
#. The sample is then removed.
#. A replacement sample is then placed into the remaining volume (i.e. it has a higher likelihood) according to :math:`\ln X_{i} \approx - (i \pm \sqrt{i})/N`, where :math:`N` is the number of samples.
#. This results in a volume contraction, which focuses on areas of higher likelihood.

The final step is to average the remaining likelihoods and to multiply it by the remaining volume varaible to get the last contribution to :math:`Z`.

The posterior weights for the :math:`i^\mathrm{th}` shell can then be written as

.. math::

   P_i = \frac{L_i[X_{i+1} - X_{i}]}{2Z}.

A density estimation method (e.g. weighted histogram) can then be used to generate the PDF.
The strength of nested sampling is that it can capture multi-modal distributions, but it can be computationally expensive.

