Fit Functions
=============

Fit functions are the mathematical representations of the data.
In quickBayes the fit functions are all classes with a very specific structure, see <insert link> for a full list of currently available options.
The first step is to create a fit function object, these are Python classes designed to make it easier for inspecting important information.
As an example we will use the gaussian function and this page will outline the various methods that you can use.
The simplest example is to just create a fit function and to then call it.

.. code-block:: python

  from quickBayes.fit_functions.gaussian import Gaussian
  import numpy as np

  function = Gaussian()
  x = np.linspace(-5, 5, 100)
  # call the function
  y = function(x, amplitude=1, x0=0, sigma=0.2)

This provides an easy way to evaluate the function.
You can also check the number of parameters by using :code:`function.N_params` to check.
However, we may want to keep the current parameters for later (e.g. they are the results of a fit).
This can be done by;

.. code-block:: python

  results = {}
  results = function.report(results, 1, 0, .2)
  print(results)

which will output

.. code-block:: python

  {'Amplitude': [1.0], 'Mean': [0], 'Sigma': [2.0]}

Repeated calls to the report method will append to the lists

.. code-block:: python

  results = function.report(results, 2, 1, .2)
  print(results)

Will output:

.. code-block:: python

  {'Amplitude': [1.0, 2], 'Mean': [0, 1], 'Sigma': [2.0, 0.2]}

It is also possible to then read the results from the report output.
This makes it easy to extract the parameters from the report (e.g. after multiple fits).
By default the method will read the first entry into the report

.. code-block:: python

  params = function.read_from_report(report)
  print(params)

Hence,

.. code-block:: python

  [1.0, 2.0, 3.0]

Or a specific index can be given

.. code-block:: python

  params = function.read_from_report(report, 1)
  print(params)

which gives

.. code-block:: python

  [2.0, 1.0, 0.2]

Similar to the `report` method is a `report_errors` method:

.. code-block:: python

   errors = [ .1, .01, .02]
   error_report = function.report_errors({}, errors, params)

The behaviour for this function is the same as :code:`report`, but some functions have non-trivial errors that will be calculated as part of this method.

When doing a fit, it can be useful to get an initial guess value :code:`guess = function.get_guess()`.
If the guess is not appropritate it can be changed by using :code:`function.set_guess([10, 4, 1])`.
Similarly the bounds for the function are given by :code:`lower, upper = function.get_bounds()`.
The bounds can be set by using :code:`function.set_bounds([1, 2, 0.4], [10, 7, 1])`.


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


Advanced Fit Functions
======================


Sometimes the problem may involve a more complicated function, such as a gaussian plus a background.
The :code:`CompositeFunction` is designed to make it easy to add multiple functions together.
To create a function that is a sum of a gaussian and a flat background, we need to use the :code:`add_function` method;

.. code-block:: python

   from quickBayes.functions.composite import CompositeFunction
   from quickBayes.functions.guassian import Gaussian
   from quickBayes.functions.BG import FlatBG

   comp = CompositeFunction
   function_1 = Gaussian()
   function_2 = FlatBG()

   comp.add_function(function_1)
   comp.add_function(function_2)


The methods report, read_from_report, report_errors will all work for the new function.
The :code:`self._prefix` mentioned earlier is how the fit functions know which function we are interested in.
This is particularly useful if we have the same function multiple times (e.g. multiple peaks).

Another advanced fit function is the :code:`ConvolutionWithResolution`.
This function takes x and y data at initialisation to represent the resolution profile, which is normalised to have an area of one.
It will then convolve the resolution profile with the sum of the functions that have been added by using :code:`add_function`.
The resolution function has an additional method :code:`update_x_range`.
This method will take the new x range and then use a spline to calculate the corresponding y values.
The normalisation is then reapplied to ensure that the resolution function is always well behaved.
For example if we want to convolve a resolution profile, :code:`rx, ry`, with a gaussian plus a delta function

.. code-block::python

   from quickBayes.functions.gaussian import Gaussian
   from quickBayes.functions.delta import Delta
   from quickBayes.functions.convolution import (
        ConvolutionWithResolution as conv)

   rx = np.linspace(-5, 5)
   ry = np.exp(-rx*rx/.8**2)

   # set a range of interest
   c_func = conv(rx, ry, -4, 4)

   function_1 = Gaussian()
   function_2 = Delta()

   c_func.add_function(function_1)
   c_func.add_function(function_2)

   # change the x range for resolution function
   new_x = np.linspace(-4., .4, 200)
   c_func.update_x_range(new_x)

The :code:`Delta` function is only well defined when used with the resolution function.

The final advanced fitting function is the :code:`QEFunction` (quasielastic function).
This has a few assumptions:

- The added peak is always the same (e.g. Lorentzian or stretch exp)
- The second parameter for the added function is always the peak centre
- The added peaks will be convolved with the resolution function
- The number of peaks in a report is less than or equal to the number of peaks in a fit function
- The elastic peak can be represented by a delta function

This function can update the resolution profile in the same way as the :code:`ConvolutionWithResolution` function.
The :code:`QEFunction` has some extra properties

- :code:`N_params` for the number of parameters in the function (including any delta function).
- :code:`N_peaks` the number of peaks that have been added to the function.

The function also has some extra methods such as :code:`add_single_function`, which defines how a new peak is added to the function.
There are also methods for setting and getting the guesses and bounds for the different parts of the function (background, delta and the function at a given index from within the sum in the convolution).
This function is an abstract class will need to be implemented, for example :code:`QlDataFunction, QSEFunction, QSEFixFunction`.
These all implement a method for adding a peak, which calls the :code:`add_single_function` method of :code:`QEFunction`.

