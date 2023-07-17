import multiprocessing
from joblib import Parallel, delayed


def parallel(max_: int, function, N: int = multiprocessing.cpu_count()):
    return Parallel(n_jobs=N)(delayed(function)(j) for j in range(max_))
