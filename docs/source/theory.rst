.. _theory:

Theory
======

In this section the key equation for quickBayes is derived in detail.
The first part is the origianl derivation as set out by Sivia.
The second part is an extension of the derivation to include distinct lines.

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




Original Derivation
-------------------

Define a model, :math:`M` as as sum of indistinguishable functions/lines

.. math::

    M = \sum_j A_j f(x, \theta)

where, :math:`A_j` is amplitude, :math:`\theta` are parameters.

The model posterior is

.. math::

   P(M|D) = P(D|M) \frac{(M)}{P(D)}

where :math:`D` is the data and :math:`P` is the probability

Assume a uniform prior :math:`P(M) = 1/N` for :math:`N` lines and :math:`P(D)` is the same for all models

.. math::

    P(M|D) \propto P(D|M)

.. math::
    :label: P(D|M)

    P(D|M) \propto \int \mathrm{d}^N\theta P(D | \theta, M) P(\theta | M)

Assume that we have a known space to investigate

.. math::
   :label: x

   x_\mathrm{min} \le x \le x_\mathrm{max}


.. math::
   :label: A

   0 \le A_j \le A_\mathrm{max} \label{A}

So the normalization prior must be the volume of the hyper cube from :math:numref:`x` and :math:numref:`A` top hat function due to asssume flat prior -> volume

.. math::
   :label: P(theta|M)

    P(\theta | M) = [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}

substituting :math:numref:`P(theta|M)` into :math:numref:`P(D|M)`

.. math::

   P(D|M) \propto [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}\int \mathrm{d}^N\theta P(D|\theta, M)

Assume that the data is subject to independent additive gaussian noise

.. math::

   P(D|\theta, M) \propto \exp\left(-\frac{\chi^2}{2}\right)

where :math:`\chi^2` is the chi squared value and is a function of :math:`\theta`

.. math::
   :label: almost

    P(D|M) \propto [(x_\mathrm{max} – x_\mathrm{min}) A_\mathrm{max}]^{-N}\int \mathrm{d}^N\theta \exp\left(-\frac{\chi^2}{2}\right)


Assume a best fit exists with corresponding best fit parameters :math:`\theta_0` and a chi squared value of :math:`\chi_\mathrm{min}^2`
Use a Taylor expansion

.. math::

   \chi^2 \approx \chi^2_\mathrm{min} + \frac{1}{2}[\underline{\theta} - \underline{\theta_0}]^\mathrm{T} \underline\nabla\ \underline\nabla \chi^2(\underline{\theta_0})[\underline{\theta} - \underline{\theta_0}]

Hence

.. math::
   :label: Taylor

   \int \mathrm{d}^N\theta \exp\left(-\frac{\chi^2}{2}\right) \approx \exp\left(-\frac{\chi^2_\mathrm{min}}{2}\right) \frac{(4\pi)^N}{\sqrt{(\mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2)) }}

where :math:`\mathrm{det}(H) = \mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2))` is the determinant of the Hessian matrix :math:`H`.
Substituting :math:numref:`Taylor` into :math:numref:`almost` and for indistinguishable lines there are :math:`N` factorial possibilities

.. math::
   :label: sivia

   P(D|M) \propto P(M|D) \propto \frac{N! (4\pi)^N }{[(x_\mathrm{max} - x_\mathrm{min})A_\mathrm{max}]^N \sqrt{\mathrm{det}(H)}} \exp\left(-\frac{\chi^2_0}{2}\right)

Taking the logs and rearranging gives

.. math::
   :label: logs

   \log{[P(D|M)]} \propto \sum_{j=1}^{N}\log{(j)} +
   N\log{(4\pi)} - N\log{([x_\mathrm{max} - x_\mathrm{min}]A_\mathrm{max})} -
   \log{(\sqrt{\mathrm{det}(H)})}  -
   \frac{\chi^2_0}{2}.

To make the equation an equality would require the addition of the normalisation for the probability, but this would be the same when comparing models with the same data set.
Hence, the term can be neglected.
The larger the probability, the larger its log value.
Hence, the better model has a larger value.
The first two terms in :math:numref:`logs` correspond to a benefit to having complexity.
This is because as the number of parameters increases, it becomes easier to fit the model to the data.
The third term is related to the prior.
If the prior is large, then little is known about the expected result.
This is then penalised as the uncertainty in the model being correct becomes larger.
The fourth term is a bit more complex as it involves the Hessian matrix.
Lets consider the best case scenario of :math:`\mathrm{det}(H) = 1`, which corresponds to a perfect fit.
Then the contribution to :math:numref:`logs` is zero.
If :math:`\mathrm{det}(H) < 1` then at least one of the eigenvalues is very small and is an indication of the model being overparameterised.
This invalidates the assumption of being at a local minima, and is accounted for in the code.
As :math:`\mathrm{det}(H)` gets larger the less likely the model.
The final term is a penality for having a poor fit, as the :math:`\chi^2` grows the likelihood decreases.


Including Unique Lines
----------------------

Define a model, :math:`M` as as sum of indistinguishable functions/lines and some other function :math:`g`

.. math::
   :label: big M

   M = \sum_i^k \alpha_i g_i(x, \underline{\theta}) + \sum_j^N A_j f(x, \underline{\theta})

where, :math:`\alpha_i` is the amplitude, :math:`\underline{\theta}` is a vector of parameters, $N$ is the number of indistinguishable lines and :math:`k` is the number of distinguishable lines.
The model posterior is

.. math::

   P(M|D) = P(D|M) \frac{(M)}{P(D)}

where :math:`D` is the data and :math:`P` is the probability
Assume a uniform prior  :math:`P(M) = 1/N` for :math:`N` lines and :math:`P(D)` is the same for all models

.. math::

   P(M|D) \propto P(D|M)

The probabilities can be split into two parts corresponding to the two terms in :math:numref:`big M`

.. math::
   P(D|M) = P(D|G + F)

where :math:`G = \sum_j \alpha_j g_j(x, \underline{\theta})` and :math:`F = \sum_j A_j f(x, \underline{\theta})`.

.. math::
    :label: P(D|G + F)

    P(D|M) \propto \int \mathrm{d}\underline{\theta} P(D | \underline{\theta}, G + F) P(\underline{\theta} | G + F)

assume that we have a known space to investigate

.. math::
   :label: x2

   x_\mathrm{min} \le x \le x_\mathrm{max}

For the :math:`F` terms:

.. math::
   :label: A2

   A_\mathrm{min} \le A_j \le A_\mathrm{max}

For the :math:`G` terms:

.. math::
   :label: alpha

   \alpha_{i_\mathrm{min}} \le \alpha_i \le \alpha_{i_\mathrm{max}}

So the normalization prior must be the volume of the hyper cube from :math:numref:`x2`, :math:numref:`A2` and :math:numref:`alpha`

.. math::
   :label: P(theta|M2)

    P(\underline{\theta} | G + F) = [(x_\mathrm{max} – x_\mathrm{min}) (A_\mathrm{max}-A_\mathrm{max})]^{-N}(x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})]^{-1}

The first part of this is just a more general version of :math:numref:`P(theta|M)`, so let :math:`\beta =  [(x_\mathrm{max} – x_\mathrm{min}) (A_\mathrm{max}-A_\mathrm{max})]^{-N}` then :math:numref:`P(theta|M2)` becomes

.. math::
   :label: P(theta|M2)2

   P(\underline{\theta} | G + F) = \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})]^{-1}


substituting :math:numref:`P(theta|M2)2` into :math:numref:`P(D|G + F)`

.. math::

   P(D|G + F) \propto \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})^{-1} \int \mathrm{d}\underline{\theta} P(D | \underline{\theta}, G + F)

Assume that the data is subject to independent additive gaussian noise

.. math::

   P(D|\underline{\theta}, G + F) \propto \exp\left(-\frac{\chi^2}{2}\right)

where :math:`\chi^2` is the chi squared value and is a function of :math:`\underline{\theta}`

.. math::
   :label: almost2

   P(D|G + F) \propto  \beta (x_\mathrm{max} – x_\mathrm{min})^{-k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})^{-1} \int \mathrm{d}\underline{\theta} \exp\left( - \frac{\chi^2}{2}\right)

Assume a best fit exists with corresponding best fit parameters :math:`\underline{\theta_0}` and a chi squared value of :math:`\chi_\mathrm{min}^2`

Use a Taylor expansion

.. math::

    \chi^2 \approx \chi^2_\mathrm{min} + \frac{1}{2}[\underline{\theta} - \underline{\theta_0}]^\mathrm{T} \underline\nabla\ \underline\nabla \chi^2(\underline{\theta_0})[\underline{\theta} - \underline{\theta_0}]

Hence

.. math::
   :label: Taylor2

   \int \mathrm{d}\underline{\theta} \exp\left(-\frac{\chi^2}{2}\right) \approx \exp\left(-\frac{\chi^2_\mathrm{min}}{2}\right) \frac{(4\pi)^{N+k}}{\sqrt{(\mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2)) }}

where :math:`\mathrm{det}(H) = \mathrm{det}(\underline{\nabla} \ \underline{\nabla} \chi^2))` is the determinant of the Hessian matrix :math:`H`.
Substituting :math:numref:`Taylor2` into :math:numref:`almost2` and for indistinguishable lines there are :math:`N` factorial possibilities

.. math::
   :label: me

   P(D|M) \propto P(M|D) \propto \frac{N! (4\pi)^{N+k}\beta }{\sqrt{H}(x_\mathrm{max} – x_\mathrm{min})^{k}\prod_i^k (\alpha_{i_\mathrm{max}}-\alpha_{i_\mathrm{max}})} \exp\left(-\frac{\chi^2_0}{2}\right)

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
   \frac{\chi^2_0}{2}
   \end{eqnarray}

If the :math:`k` distinguishable lines are the same for all models being considered, then the :math:`k\log{(x_\mathrm{max} - x_\mathrm{min})}`, :math:`k\log{(4\pi)}` and :math:`\sum_i^k
\log{(\alpha_{i_\mathrm{max}}-
\alpha_{i_\mathrm{max}})}` terms can be neglected as they just add a constant offset. Hence,

.. math::

   \log{[P(D|M)]} \propto \sum_{j=1}^{N}\log{(j)} +
   N\log{(4\pi)} + \log{(\beta)} -
   \log{(\sqrt{H})}  -
   \frac{\chi^2_0}{2}

In the case of positive definite amplitudes :math:`A_\mathrm{min} = 0` and substituting in for :math:`\beta` this reduces to :math:numref:`logs`.
Alternatively, substituting :math:numref:`me` into the odds ratio would lead to the terms corresponding to the distinguishable lines cancelling out.
So they can be neglected, this might happen in the case of a linear background term for all of the models.
