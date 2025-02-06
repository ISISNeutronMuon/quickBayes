.. _direct:


Comparison with direct Bayesian methods
---------------------------------------

In this section quickBayes will be compared to methods that sample the Bayesian posterior directly (:ref:`MCMC <MCMC>` and :ref:`nested sampling <nest>`).
These methods that sample the posterior distribution directly are computationally expensive and normally slow.
From these posterior distributions different models can be compared, using the odds factor.
The `quickBayes` package is significantly less computationally demanding than its direct Bayesian counterparts.
However, this is at the cost of only calculating the likelihood for comparing different models.
The posterior distributions are never calculated explicitly.
Hence, why its not as computationally intensive.

For both the Bayesian methods and `quickBayes` the equation of interest is the probability of the data given the model.
Therefore, `quickBayes` is attempting to solve the exact same problem as other Bayesian methods by making simplifing assumptions.

The `quickBayes` package is best used when the user just wants to know which model is most likely.
If the user wants to know the posterior PDFs then one of the Bayesian methods would be more appropriate.

