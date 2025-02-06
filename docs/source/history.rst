.. _history:

History of quickBayes
=====================

The quickBayes packages started as a replacement for the quasielasticbayes package, which was a Fortran code base that was exposed to Python.
In the development of quickBayes the code became more modular and flexible.
This allows quickBayes to be applied to any problem that requires model selection.

To verify the results of quickBayes it was compared to the quasielasticbayes calculations.
For lorentzian fits of QENS data it produced results that were close enough to quasielasticbayes.
However, for a stretched exponential quickBayes gave different results for the Full Width Half Max (FWHM)
This is shown by the figure below.
The quickBayes results agree with the FWHM values for fitting a single Lorentzian ('QL'), suggesting that the results are correct.
Whereas the quasielasticbayes package has a divergence in the FWHM for low :math:`Q` values.

.. figure:: /images/qse_cf.png
   :alt: qse_cf.png
   :width: 400px
   :align: center
