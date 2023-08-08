Parallel
========

The parallel function is included to make it easier to do similar calculations faster.
The parallel function uses `Joblib <https://joblib.readthedocs.io/en/stable/index.html>`_.

A simple example is:

.. code-block:: python

  import time
  import numpy as np
  from quickBayes.utils.parallel import parallel


  def do_thing(j):
          time.sleep(0.5)
          return -j

  data = np.zeros(100)

  start = time.time()
  data = Parallel(list(range(100)), do_thing)
  print(data)
  print(time.time() - start)

Notice that the output is a list of the results from the :code:`do_thing` function and that the :code:`range` has to be converted to a list.

If you have multiple inputs then you can pass them by using tuples.

.. code-block:: python

   from quickBayes.utils.parallel import parallel

   def function(input):
       i = input[0]
       j = input[1]
       return i*j

   inputs = [ (k, -k) for k in range(10)]
   result = parallel(inputs, function)
   print(result)

There will be some occasions when some of the parameters are the same for all of the parallel functions.
For this case the :code:`partial` function of :code:`functools` can be used:

.. code-block:: python

   from functools import partial
   from quickBayes.utils.parallel import parallel

   def function(input, offset):
       i = input[0]
       j = input[1]
       return i*j - offset

   inputs = [ (k, -k) for k in range(10)]
   const = 3.
   partial_function = partial(function, offset=const)
   result = parallel(inputs, partial_function)
   print(result)

The :code:`parallel` function can even be used with a method from a class.

.. code-block:: python

   from functools import partial
   from quickBayes.utils.parallel import parallel

   class test_parallel(object):
       def __init__(self, offset):
           self._offset = offset

       def function(self, input, constant):
           i = input[0]
           j = input[1]
           return i*j - constant*self._offset

       def parallel_function(self, values, constant):
           partial_function = partial(self.function, constant=constant)
           return parallel(values, partial_function)

   inputs = [ (k, -k) for k in range(10)]
   const = 3.
   test = test_parallel(2.)
   result = test.parallel_function(inputs, const)
   print(result)

It is important to note that the method call is to a wrapper of the method function we want to run in parallel.

