Loglikelihood
=============

Calculate the unnormalised logliklihood.
log_10 probabilities -> priors
The equation for the probability is taken from "Data Analysis A bayesian tutorial second edition", by D. S. Sivia equation 4.20 page 88:

.. math::
    P(M|D_k,I) \propto \frac{N! (4\pi)^N}{\beta^N\sqrt{\det(H)}}\exp{(-\frac{\chi^2}{2})},

where :math:`\\beta = (x_{max}-x_{min})A_{max}`, :math:`H` is the Hessian matrixWe want this as a logliklihood, so take :math:`\\log_{10}` and use:

.. math::
    \log_{10}(\exp{\alpha}) = \ln(\exp{\alpha})\log_{10}(\exp{1}) \\
    \log_{10}(\exp({\alpha}) = \alpha\log_{10}(\exp{1}) \\
    N! = \prod_{i=0}^{N} i => \log{N!} = \sum_{i=0}^N \log(i) \\
    \log(ab) = \log(a) + \log(b) \\
    \log(\frac{a}{b}) = \log(a) - \log(b) \\
    \log(a^N) = N\log(a)

to get:

.. math::
    \sum_{j=1}{N}\log(j) + N\log(4\pi) - 0.5\chi^2\log(\exp(1))
    - N\log(\beta) - 0.5\log(\det(H))



