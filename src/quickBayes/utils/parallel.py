import multiprocessing
from joblib import Parallel, delayed
from collections.abc import Callable


def parallel(items: list, function: Callable,
             N: int = multiprocessing.cpu_count()):
    """
    This is a wrapper of the joblib Parallel function.
    It will run the function over multiple cores and then return the result.
    Note that the function must take the looped value as its only input.
    Use threads as the default does not work with Mantid.
    :input items: the list to loop over
    :input function: the function to run in parallel
    :input N: the number of process to use
    :return a list of the outputs from function. If multuple outputs
    from function then the first index is for the loop value and the
    second index is for the item from function.
    """
    return Parallel(n_jobs=N,
                    prefer="threads")(delayed(function)(j) for j in items)
