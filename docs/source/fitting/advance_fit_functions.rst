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

.. code-block:: python

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

