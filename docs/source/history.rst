History of quickBayes
=====================

The quickBayes packages started as a replacement for the quasielasticbayes package, which was a Fortran code base that was exposed to Python.
In the development of quickBayes the code became more modular and flexable. 
This allows the code to be applied to any probelm that requires model selection. 

To verify the results of quickBayes it was compared to the quasielasticbayes calculations. 
For lorentzian fits it produced the close enough results.
However, for a stretched exponential it gave different results for the FWHM.
The stretched exponential results for the FWHM are different to quasielasticbayes, as shown by the figure below.
However, the new results agree with the FWHM values for fitting a single Lorentzian ('QL').
The quickBayes method provides FWHM results that are comparable for all :math:`Q` values (green and black data), unlike the original code that has a divergence for low :math:`Q` value.

.. figure:: /images/qse_cf.png
   :alt: qse_cf.png
   :width: 400px
   :align: center
