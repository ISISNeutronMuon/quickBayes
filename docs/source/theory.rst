.. _theory:

The theory behind quickBayes
============================

In this section the key equation for quickBayes is derived in detail.
The quickBayes method makes a series of assumptions to reduce :ref:`the full PDF evaluation <eq_int>` to a single analytic equation.
This will follow the derivation of `Sivia et al <https://www.sciencedirect.com/science/article/pii/092145269290036R?via=ihub>`_.
The key assumptions are:

- The model can be written as a series of indistinguishable lines (i.e. the same repeated function)
- The lines can be written as an amplitude multiplied by some function
- The prior probabilities are flat across the domain of interest
- The normalisation is just the value of the hyper volume for the parameters
- The :math:`P(D|\underline\theta, M)` can be represented as gaussians and approximated by a first order Taylor expansion

In the second part the original derivation has been extended to remove the first assumption, allowing for unique lines (functions) to be included.

When used by the quickBayes package the chi squared values are calculated using either `scipy <https://scipy.org/>`_ or `gofit <https://ralna.github.io/GOFit/_build/html/index.html>`_.
However, it is possible to add new fitting engines if these are not sufficient.

Original Derivation
-------------------

We want to represent the data with a model, but there is some uncertainty about which model to use.
We will assume that the models are similar, with the only difference being the number of functions being used to describe the data.
Hence, the model, :math:`M`, can be written as as sum of indistinguishable functions/lines

.. math::

    M = \sum_j A_j f(x, \theta)

where, :math:`A_j` is amplitude, :math:`\theta` are parameters and :math:`f` is the function.

From :math:numref:`bayes` we know that the model posterior can be written as

.. math::

   P(M|D) = P(D|M) \frac{P(M)}{P(D)}

where :math:`D` is the data and :math:`P` is the probability.
If we assume a uniform prior :math:`P(M) = 1/N` for :math:`N` lines and that :math:`P(D)` is the same for all models then

.. math::
   :label: prior

    P(M|D) \propto P(D|M).

From :math:numref:`eq_int` we know that :math:`P(D|M)` can be written as an intergral to get,

.. math::
    :label: P(D|M)

    P(D|M) \propto \int \mathrm{d}^N\theta P(D | \theta, M) P(\theta | M).

To make the above equation tractable, let's assume that we know that the probability is zero outside of some known space such that

.. math::
   :label: x

   x_\mathrm{min} \le x \le x_\mathrm{max}

and

.. math::
   :label: A

   0 \le A_j \le A_\mathrm{max}.

This can just be thought of as defining the prior.
With no additional information on the prior, we assume that it is just flat within this space.
Therefore, the prior probability is just the volume of the hyper cube from :math:numref:`x` and :math:numref:`A`

.. math::
   :label: P(theta|M)

    P(\theta | M) = [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}.

Substituting :math:numref:`P(theta|M)` into :math:numref:`P(D|M)` yields

.. math::
   :label: P(DM)

   P(D|M) \propto [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}\int \mathrm{d}^N\theta P(D|\theta, M).

To continue simplifing we will assume that the data is subject to independent additive gaussian noise.
Hence,

.. math::
   :label: P(D|theta,M)_exp

   P(D|\theta, M) \propto \exp\left(-\frac{\chi^2}{2}\right)

where :math:`\chi^2` is the chi squared value and is a function of the fit parameters :math:`\theta`.
Substituting this into :math:numref:`P(DM)` gives

.. math::
   :label: almost

    P(D|M) \propto [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}\int \mathrm{d}^N\theta \exp\left(-\frac{\chi^2}{2}\right).

The next step is to assume that a best fit exists, and that the corresponding best fit parameters are :math:`\theta_0` and the chi squared value is :math:`\chi_\mathrm{min}^2`.
A Taylor expansion of the chi squared yields

.. math::
   :label: chi2_expansion

   \chi^2 \approx \chi^2_\mathrm{min} + \frac{1}{2}[\underline{\theta} - \underline{\theta_0}]^\mathrm{T} \underline\nabla\ \underline\nabla \chi^2(\underline{\theta_0})[\underline{\theta} - \underline{\theta_0}]

and the integral can then be written as

.. math::
   :label: Taylor

   \int \mathrm{d}^N\theta \exp\left(-\frac{\chi^2}{2}\right) \approx \exp\left(-\frac{\chi^2_\mathrm{min}}{2}\right) \frac{(4\pi)^N}{\sqrt{(\mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2)) }},

where :math:`\mathrm{det}(H) = \mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2))` is the determinant of the Hessian matrix :math:`H`.
Substituting :math:numref:`Taylor` into :math:numref:`almost` and for :math:`N` indistinguishable lines there are :math:`N!` possibilities

.. math::
   :label: sivia

   P(D|M) \propto P(M|D) \propto \frac{N! (4\pi)^N }{[(x_\mathrm{max} - x_\mathrm{min})A_\mathrm{max}]^N \sqrt{\mathrm{det}(H)}} \exp\left(-\frac{\chi^2_0}{2}\right).

Finally, by taking the logs and rearranging this equation gives

.. math::
   :label: logs

   \log{[P(D|M)]} \propto \sum_{j=1}^{N}\log{(j)} +
   N\log{(4\pi)} - N\log{([x_\mathrm{max} - x_\mathrm{min}]A_\mathrm{max})} -
   \log{(\sqrt{\mathrm{det}(H)})}  -
   \frac{\chi^2_0}{2}.

To make the equation an equality would require the addition of the normalisation for the probability, but this would be the same when comparing models with the same data set.
Hence, the term can be neglected.
As the probability increases then the log will become more positive.
Hence, the better model has a larger value for :math:numref:`logs`.
The first two terms in :math:numref:`logs` correspond to a benefit to having complexity.
This is because as the number of parameters increases, it becomes easier to fit the model to the data.
The third term is related to the prior.
If the prior is large, then little is known about the expected result.
This is then penalised as the uncertainty in the model being correct becomes larger.
The fourth term is a bit more complex as it involves the Hessian matrix.
Lets consider the best case scenario of :math:`\mathrm{det}(H) = 1`, which corresponds to a perfectly behaved model.
Then the contribution to :math:numref:`logs` is zero.
If :math:`\mathrm{det}(H) < 1` then at least one of the eigenvalues is very small and is an indication of the model being overparameterised.
This invalidates the assumption of being at a local minimum, and quickBayes will automatically add an additional penality if this occurs.
As :math:`\mathrm{det}(H)` gets larger the less likely the model is to be correct.
The final term is a penality for having a poor fit, as the quality of the fit decreases so does the likelihood of the model.

Including Unique Lines
----------------------

Sometimes we will want to determine the best model when there are distinguishable lines within the model.
For example, selecting if the background is flat or linear.
In this section we will show how to handle this sort of model selection analytically by using a series of approximations and assumptions.
Since the derivation is very similar to the above, just the key changes will be highlighted here.

Lets define a model, :math:`M` as a sum of indistinguishable functions/lines and some other functions :math:`g_i`

.. math::
   :label: big M

   M = \sum_i^k \alpha_i g_i(x, \underline{\theta}) + \sum_j^N A_j f(x, \underline{\theta})

where, :math:`\alpha_i` is the amplitude of the :math:`i^\mathrm{th}` distinguishable function, :math:`\underline{\theta}` is a vector of parameters, :math:`N` is the number of indistinguishable lines and :math:`k` is the number of distinguishable lines.
Once again the model posterior is

.. math::

   P(M|D) = P(D|M) \frac{(M)}{P(D)}.

Assuming that the prior is uniform yields :math:numref:`prior`, but the evidence can then be split into two parts corresponding to the two terms in :math:numref:`big M`

.. math::

   P(D|M) = P(D|G + F),

where :math:`G = \sum_j \alpha_j g_j(x, \underline{\theta})` and :math:`F = \sum_j A_j f(x, \underline{\theta})`.
Hence, :math:numref:`P(D|M)` can be written as

.. math::
    :label: P(D|G + F)

    P(D|M) \propto \int \mathrm{d}\underline{\theta} P(D | \underline{\theta}, G + F) P(\underline{\theta} | G + F).

We then assume that the bounds for the prior are known, with the :math:`x` values being

.. math::
   :label: x2

   x_\mathrm{min} \le x \le x_\mathrm{max}

and the amplitudes of the :math:`F` terms are

.. math::
   :label: A2

   A_\mathrm{min} \le A_j \le A_\mathrm{max}.

For the distinguishable lines (:math:`G` terms) the bounds for the :math:`i^\mathrm{th}` term can be written as

.. math::
   :label: alpha

   \alpha_{i_\mathrm{min}} \le \alpha_i \le \alpha_{i_\mathrm{max}}.

The prior is still given by the volume of the hyper cube from :math:numref:`x2`, :math:numref:`A2` and :math:numref:`alpha`, which gives

.. math::
   :label: P(theta|M2)

    P(\underline{\theta} | G + F) = [(x_\mathrm{max} – x_\mathrm{min}) (A_\mathrm{max}-A_\mathrm{max})]^{-N}(x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})]^{-1}.

The first part of this is just :math:numref:`P(theta|M)`.
To simplify the notation let :math:`\beta =  [(x_\mathrm{max} – x_\mathrm{min}) (A_\mathrm{max}-A_\mathrm{max})]^{-N}`, which is the contribution to the prior for the distinguishable lines, then :math:numref:`P(theta|M2)` becomes

.. math::
   :label: P(theta|M2)2

   P(\underline{\theta} | G + F) = \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})]^{-1}.


Substituting :math:numref:`P(theta|M2)2` into :math:numref:`P(D|G + F)` gives

.. math::

   P(D|G + F) \propto \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})^{-1} \int \mathrm{d}\underline{\theta} P(D | \underline{\theta}, G + F).

Once again we can assume that the data is subject to independent additive gaussian noise

.. math::

   P(D|\underline{\theta}, G + F) \propto \exp\left(-\frac{\chi^2}{2}\right).

Hence,

.. math::
   :label: almost2

   P(D|G + F) \propto  \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})^{-1} \int \mathrm{d}\underline{\theta} \exp\left( - \frac{\chi^2}{2}\right)

and we can assume that a best fit exists with corresponding best fit parameters :math:`\underline{\theta_0}` and a chi squared value of :math:`\chi_\mathrm{min}^2`.
The Taylor expansion in :math:numref:`Taylor` can then be written as

.. math::
   :label: Taylor2

   \int \mathrm{d}\underline{\theta} \exp\left(-\frac{\chi^2}{2}\right) \approx \exp\left(-\frac{\chi^2_\mathrm{min}}{2}\right) \frac{(4\pi)^{N+k}}{\sqrt{(\mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2)) }}.

Substituting :math:numref:`Taylor2` into :math:numref:`almost2` and including a factor of :math:`N!` for the possibilities of :math:`N` indistinguishable lines

.. math::
   :label: me

   P(D|M) \propto P(M|D) \propto \frac{N! (4\pi)^{N+k}\beta }{\sqrt{H}(x_\mathrm{max} – x_\mathrm{min})^{k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})} \exp\left(-\frac{\chi^2_0}{2}\right).

Taking the log of this expression and rearranging yields

.. math::
   :nowrap:

   \begin{eqnarray}
   \log{[P(D|M)]} \propto \sum_{j=1}^{N}\log{(j)} +
   (N+k)\log{(4\pi)} + \log{(\beta)} -
   \log{(\sqrt{H})} \\ -
   k\log{(x_\mathrm{max} - x_\mathrm{min})}
   - \sum_i^k
   \log{(\alpha_{i_\mathrm{max}}-
   \alpha_{i_\mathrm{max}})} -
   \frac{\chi^2_0}{2}.
   \end{eqnarray}

If the :math:`k` distinguishable lines are the same for all models being considered, then the :math:`k\log{(x_\mathrm{max} - x_\mathrm{min})}`, :math:`k\log{(4\pi)}` and :math:`\sum_i^k
\log{(\alpha_{i_\mathrm{max}}-
\alpha_{i_\mathrm{max}})}` terms can be neglected as they just add a constant offset.
Hence, the above equation simplifies to

.. math::

   \log{[P(D|M)]} \propto \sum_{j=1}^{N}\log{(j)} +
   N\log{(4\pi)} + \log{(\beta)} -
   \log{(\sqrt{H})}  -
   \frac{\chi^2_0}{2}.

In the case of positive definite amplitudes :math:`A_\mathrm{min} \ge 0` and substituting in for :math:`\beta` this reduces to :math:numref:`logs`.
Alternatively, substituting :math:numref:`me` into the odds ratio would lead to the terms corresponding to the distinguishable lines cancelling out.
This happens when the models all include the same background term (e.g. flat) and then only differ by the number of distinguishable lines.
