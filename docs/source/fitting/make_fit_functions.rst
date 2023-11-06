Making A New Fit Function
=========================

All fit functions must inherit from :code:`BaseFitFunction`, this defines the basic methods of a fit function.
Each class will need the following methods:

- :code:`__init__` to define the setup of the function
- :code:`__call__` to define how to evaluate the function
- :code:`read_from_report` to define how to extract values from a dictionary
- :code:`report` to define how to add the parameters to a dictionary
- :code:`report_errors` to define how to add the parameter errors to a dictionary.

It is then strongly recommended that the class also inludes properties corresponding to each of the variables.
These methods will simply return the name of the variable, this prevents issues with typos.
For example

.. code-block:: python

   @property
   def amplitude(self) -> str:
       """
       :return string for amplitude name
       """
       return f"{self._prefix}Amplitude"

this example uses an f string and the :code:`self._prefix` to handle Advanced fitting functions.

