Parallel
========

The parallel function is included to make it easier to do similar calculations faster.
The parallel function uses `Joblib <https://joblib.readthedocs.io/en/stable/index.html>`_.

A simple example is:
.. code-block:: python
  :linenos:
        import time
        import numpy as np
        import multiprocessing
        from multiprocessing import Pool, freeze_support
        from multiprocessing import Pool
        from joblib import Parallel, delayed


        start = time.time()

        def do_thing(j):
                time.sleep(0.5)
                return -j

        data = np.zeros(100)

        data = Parallel(list(range(100)), do_thing)
        print(data)
        print(time.time() - start)

Notice that the output is a list of the results from the `do_thing` function.
