.. _AIC:

Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC)
---------------------------------------------------------------------------

Both the `AIC <https://en.wikipedia.org/wiki/Akaike_information_criterion#>`_ and `BIC <https://en.wikipedia.org/wiki/Bayesian_information_criterion>`_ are methods for determining which model best suits the data.
These methods have their origin in information theory.
The key idea is to reduce the amount of information lost by the model when describing the data.

Lets assume that a function, :math:`g`, exists that perfectly describes the observed data.
The Kullback-Leibler distance then defines the information as

.. math::

   I = E\left [\ln\left \{\frac{g(x)}{P(x|\underline{\theta})}\right\}\right],

where :math:`x` are the observed data points, and :math:`E` is a functional defined as

.. math::

   E[y(x)] = \int g(x) y(x) \mathrm{dx}.

However, the exact form of the functional does not impact the derivation.
As a result the information can be written as

.. math::

   I = E\left[ \ln{\{g(x)\}} - \ln\{f(x | \underline{\theta})\}\right].

It is clear that the way to minimize the information loss is to minimise the argument for the functional,

.. math::
   :label: AIC_derivation

   a = \ln{\{g(x)\}} - \ln\{P(x | \underline{\theta}\}.

Hence, the best model is the one with the lowest value, this reamins true for the AIC and BIC.
The only part that depends on the model in equation :math:numref:`AIC_derivation` is :math:`P(x|\underline{\theta})`, which is liklihood of the parameters (and model) given the data and will be denoted by :math:`\mathcal{L}(\underline{\theta}|x)`.
This allows equation :math:numref:`AIC_derivation` to be written as

.. math::
   :label: AIC_no_approx

   a = \ln{\{g(x)\}} - \ln{\{\mathcal{L}(\underline{\theta}|x)\}}.

When comparing models to the same data set, the first term will be a constant and only the second term depends on the choice of model.
However, the exact form of :math:`g(x)` is not known and as a result it must be estimated.
The AIC and BIC use different approximations for this first term.
Both the AIC and BIC multiply equation :math:numref:`AIC_no_approx` by a factor two (for historic reasons).
To derive the AIC some `statistical arguments <https://ieeexplore.ieee.org/document/1100705>`_ are made that :math:`2\ln{\{g(x)\}} = 2k` to get

.. math::

   \mathrm{AIC} = 2k - 2\ln{\{\mathcal{L}(\underline{\theta}|x)\}},

where :math:`k` are the number of parameters in the model.
This approximation only holds true for large sample sizes, which results in it giving preference to overparameterised models for a small number of data points.
A more sophisticated version of the AIC has been developed to account for small sample sizes

.. math::

   \mathrm{AIC_c} = \mathrm{AIC} + \frac{2k^2 + 2k}{n - k - 1},

where :math:`n` is the number of data points.
It is clear that in the limit of infinite data this just reduces to the standard AIC.
Whereas the BIC `shows that <https://www.jstor.org/stable/2958889>`_ (via a different derivation) :math:`2\ln{\{g(x)\}} = k\ln{(n)}`, which gives

.. math::

   \mathrm{BIC} = k\ln{n} - 2\ln{\{\mathcal{L}(\underline{\theta}|x)\}}.


To better understand these methods it is worth considering the case of a gaussian likelihood function

.. math::

   \mathcal{L}(\underline{\theta}|x) = \frac{1}{\sigma\sqrt{2\pi}}\exp\left(-\frac{\sum_j (y_j - h(\underline{\theta}, x))^2}{2\sigma^2}\right),

where :math:`h` is the model (fitting function) being used to describe the data, :math:`y_j` is the observed :math:`j^\mathrm{th}` data point and :math:`\sigma` is the uncertainty.
This means that the likelihood can be written as

.. math::

   \mathcal{L}(\underline{\theta}|x) = C\exp\left(-\frac{\chi^2}{2}\right),

where :math:`\chi^2` is the chi squared value from linear least squares and :math:`C` is a constant term.
Since only differences are important, equation :math:numref:`AIC_no_approx` can be written as

.. math::

   2a = 2\ln{\{g(x)\}} - \chi^2.

For both the AIC and BIC the first term is the same if the comparing two models with the same number of parameters against the same data set.
Hence, the best AIC and BIC is just the model with the lowest :math:`\chi^2` value.


Comparison with AIC and BIC
---------------------------

The AIC and BIC are both (relatively) simple equations for calculating the most likely model.
This is similar to the ethos behind the `quickBayes` package.
However, the AIC and BIC both originate from information theory, while `quickBayes` starts from the probability of the data given the model.
To explore this distinction we will consider a pair of models

.. math::
   :label: cf_f_def

   M_N(x, \underline{\theta}) = \sum_{j}^N f(x, \underline{\theta}),

where the repeated function :math:`f` is repeated :math:`N` times, with the parameters :math:`\underline{\theta}`.
When increasing the number of lines by one, the number of fitting parameters will increase by :math:`k`.
To compare two AIC's we can subtract two neighbouring models from each other

.. math::

   \Delta \mathrm{AIC} = \mathrm{AIC}_{N+1} - \mathrm{AIC}_N,

where the :math:`\mathrm{AIC}_N` is an AIC with :math:`N` functions.
Assuming a gaussian distribution, this can be simplified to

.. math::
   :label: Delta_AIC

   \Delta AIC = 2k + \chi_N^2 - \chi_{N+1}^2,

where :math:`\chi_m^2` is the chi squared value for a model with :math:`N` functions.
Similarly, the change in BIC due to two neighbouring models can be written as

.. math::

   \Delta \mathrm{BIC} = \mathrm{BIC}_{N+1} - \mathrm{BIC}_N,

where the :math:`\mathrm{BIC}_N` is an BIC with :math:`N` functions.
It can be shown that for a gaussian distribution,

.. math::
   :label: Delta_BIC

   \Delta BIC = 2k\ln{(n)} + \chi_N^2 - \chi_{N+1}^2,

where :math:`n` is the number of data points.
The interpretation of equations :math:numref:`Delta_AIC` and :math:numref:`Delta_BIC` are similar.
When comparing the AIC/BIC it is the value which is smaller that is most likely.
Hence, if equations :math:numref:`Delta_AIC` or :math:numref:`Delta_BIC` are negative then the model with :math:`N+1` functions is prefered.
Alternatively, if the value is positive then less function (i.e. :math:`N`) are prefered.
For both equations :math:numref:`Delta_AIC` and :math:numref:`Delta_BIC` they have a cost term for adding an extra function and then a difference in the chi squared values.
If the number of parameters per function (:math:`k`) is much smaller than the difference in the chi squared values, then this is equivalent to just comparing the goodness of fit.

The main equation for quickBayes (equation :math:numref:`logs`) can be written as

.. math::

   \ln{[P(D|M_N)]} = C + \sum_{j=1}^{N}\ln{(j)} +
   N\ln{(4\pi)} - N\ln{([x_\mathrm{max} - x_\mathrm{min}]A_\mathrm{max})} -
   \ln{(\sqrt{\det{H}})}  -
   \frac{\chi^2}{2},

where :math:`\chi^2` is at the minimum, :math:`M_N` is the model with :math:`N` functions and :math:`C` is a normalisation constant.
The normalisation constant will be the same for all of the models, so when taking the difference it will cancel out.
Lets define the difference between two neighbouring models to be

.. math::

   \Delta = \ln{[P(D|M_{N+1})]} - \ln{[P(D|M_N)]}.

When using :math:`\Delta` to determine the best model, some care is needed.
The probabilities should always be less than one, so the logs are negative, and the as a model becomes more probable the value gets closer to zero.
Hence, a negative result for :math:`\Delta` means that :math:`N` functions are prefered and a positive result means :math:`N+1` functions give the more likely fit.
Substituting in :math:numref:`logs` into the defintion of :math:`\Delta` yields,

.. math::
   \Delta = \ln(4\pi) - \ln{([x_\mathrm{max} - x_\mathrm{min}]A_\mathrm{max})} - \ln{(\sqrt{\det{H_{N+1}}})} - \frac{\chi^2_{N+1}}{2} + \ln{(\sqrt{\det{H_{N}}})} + \frac{\chi^2_{N}}{2},

which can be rearranged to

.. math::
   :label: Delta_qb

   \Delta = \ln(4\pi) - \ln{([x_\mathrm{max} - x_\mathrm{min}]A_\mathrm{max})} + \ln{\left(\frac{\sqrt{\det{H_{N}}}}{\sqrt{\det{H_N+1}}}\right)} + \frac{1}{2}(\chi^2_{N} - \chi^2_{N+1}).

The first term is a benefit to using more functions and the second term is the cost of the prior, as discussed preveously.
The third term is related to the Hessian matrix.
The final term is the difference between the chi squared values, similar to the differences in the AIC and BIC.

Equations :math:numref:`Delta_AIC` and :math:numref:`Delta_BIC` have a term that penalises complexity (more functions) and then a term that is the comparison between the goodness of the fits for the two models.
Whereas, equation :math:numref:`Delta_qb` instead has a term that encourages complexity and a term that penalises poor prior knowledge.
The comparison for the goodness of fit also takes into account the Hessian matrices for the models.
