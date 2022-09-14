import numpy as np
from scipy import linalg
import os
from libc.math cimport log10, round, ceil, pow, abs
import cython
from libc.math cimport sin
cimport numpy as np


np.import_array()


"""
This file contains methods and classes to make converting Fortran to Python easier.
This is achived by this class making adjustments so that the syntax of the original
code is changed by as little as possible

"""

def update_vec(np.ndarray[np.float_t] new_vec, np.ndarray[np.float_t] vec):
    vec[:len(new_vec)] = new_vec
    return vec



cpdef round_sig(float x,int sig=8, float small_value = 1.0e-9):
    if x == 0.0:
        return 0.0
    cdef float a = abs(x)
    cdef float l = log10(a)
    cdef float power = sig - ceil(log10(abs(x)))
    cdef float factor = pow(10, power)

    result = round(x*factor)
    return result/factor


cpdef get_range(start, end, dx=1):
    return range(start, int(end+1*np.sign(dx)),int(dx))


def NINT(value):
    return int(round(value))


def CMPLX(A,B):
    return np.cdouble(complex(A,B))


cpdef complex_zeros(n):
    return np.asarray([complex(0,0) for _ in range(n)])


cdef zeros(n, complex=False):
    if complex:
        return complex_zeros(n)
    else:
        return np.zeros(n)


class Vec(object):
    """
    Vec is a simple method for making a
    numpy array look like a fortran array
    - e.g. start at 1 instead of 0
    """
    def __init__(self, n,complex=False):
        self._vec = zeros(n, complex)
        self._complex = complex

    def __call__(self,j):
        return self._vec[j-1]

    def _print(self):
        print(self._vec)

    def set(self, j, value):
        self._vec[j-1] = value

    def output(self):
        return self._vec

    def output_range(self,start=1,end=2):
        return self._vec[start-1:end+1]

    def fill(self, value, n,start=1):
        self._vec[start-1:n+start-2] = value

    def get_from(self, start=1):
        return self._vec[start-1:]

    def copy(self, vector, start_i=1):
       self._vec[start_i-1:start_i-1+len(vector)] = vector

    def dimensions(self):
        return 1


class BoolVec(Vec):
    """
    A boolean vector that looks like fortran
    """
    def __init__(self, n,complex=False):
        self._vec = [True for _ in range(n)]


class Matrix_2D(Vec):
    """
    Matrix_2D is a class for a 2D matrix.
    Under the hood this is a vector to
    replicate the behaviour of fortran,
    which allows for a silent conversion
    between a matrix and a vector. It
    also allows for fortran like indexing.
    """
    def __init__(self, m,n,complex=False):
        super().__init__(n*m, complex)
        self._m = m
        self._n = n

    def _index_for_vec(self, i, j):
        return (j-1)*self._m+i-1

    def fill(self, value, n, start_i=1,start_j=1):
        k = self._index_for_vec(start_i,start_j)
        super().fill(value,n,k+1)

    def __call__(self,i,j):
        return self._vec[self._index_for_vec(i,j)]

    def set(self, i,j, value):
        self._vec[self._index_for_vec(i,j)] = value

    def _print(self):
        print("i, j, value")
        for j in range(self._n):
            for i in range(self._m):
                k = self._index_for_vec(i+1,j+1)
                print(i+1,j+1,self._vec[k])

    def output_from(self, i,j):
        k = self._index_for_vec(i,j)
        return self._vec[k:]

    def output_range(self, i,j,end):
        k = self._index_for_vec(i,j)
        return self._vec[k:k+end]

    def output_col(self, j):
        k = self._index_for_vec(1,j)
        return self._vec[k:k+self._m]

    def mat(self):
        mat = np.asarray([np.asarray(self.output_col(j)) for j in get_range(1,self._n)])
        return mat

    def copy(self, vec, start_i=1, start_j=1):
       k = self._index_for_vec(start_i,start_j)
       self._vec[k:k+len(vec)] = vec

    def output_as_vec(self):
        tmp = Vec(self._m*self._n)
        tmp.copy(self.output())
        return tmp

    def resize(self, new_m, new_n):

        max_m = self._m
        max_n = self._n
        if new_m < self._m or new_n < self._n:
            max_m = new_m
            max_n = new_n
        new_vec = zeros(new_n*new_m, self._complex)
        k2 = 0
        for j in range(max_n):
            for i in range(max_m):
                k = self._index_for_vec(i+1,j+1)
                new_vec[k2] = self._vec[k]
                k2+=1
            k2+= (new_m -self._m)
        self._m = new_m
        self._n = new_n
        self._vec = new_vec

    def dimensions(self):
        return 2


cdef class C_Vec(object):
    """
    C_Vec is a Cython version of Vec
    for speed.
    """
    cdef double[:] _vec

    def __cinit__(self,int n):
        self._vec = np.zeros(n)

    def __call__(self,j):
        return self._vec[j-1]

    cpdef set(self,size_t j,double value):
        self._vec[j-1] = value

    cpdef output(self):
        return np.array(self._vec)

    cpdef output_c(self):
        return self._vec

    cpdef output_range(self,size_t start=1,size_t end=2):
        cdef double[:] tmp = self._vec[start-1:end+1]
        return np.array(tmp)

    cpdef fill(self,float value,size_t n, size_t start=1):
        self._vec[start-1:n+start-2] = value

    cpdef get_from(self,size_t start=1):
        cdef double[:] tmp = self._vec[start-1:]
        return np.array(tmp)

    cpdef copy(self,np.ndarray[double] vector, size_t start_i=1):
       cdef int k
       for k in range(len(vector)):
            self._vec[start_i+k-1] = vector[k]

    cpdef dimensions(self):
        return 1


cdef class C_Matrix_2D(object):
    """
    C_Matrix_2D is a Cython version
    of Matrix_2D for speed
    """
    cdef int _m, _n
    cdef double[:] _vec

    def __init__(self,int m, int n):
        self._vec = np.zeros(n*m)
        self._m = m
        self._n = n

    cdef _index_for_vec(self,size_t i,size_t j):
        return (j-1)*self._m+i-1

    cpdef fill(self,float value,int n,size_t start_i=1, size_t start_j=1):
        cdef size_t k = self._index_for_vec(start_i,start_j)
        self._vec[k:n+k-1] = value

    def __call__(self,size_t i, size_t j):
        cdef size_t k =self._index_for_vec(i,j)
        return self._vec[k]

    cpdef set(self, i,j, value):
        cdef size_t k =self._index_for_vec(i,j)
        self._vec[k] = value

    cpdef output_from(self,size_t i, size_t j):
        cdef size_t k = self._index_for_vec(i,j)
        cdef double[:] tmp = self._vec[k:]
        return np.array(tmp)

    cpdef output_range(self,size_t i,size_t j,size_t end):
        cdef size_t k = self._index_for_vec(i,j)
        cdef double[:] tmp = self._vec[k:k+end]
        return np.array(tmp)

    cpdef output(self):
        return np.array(self._vec)

    cpdef output_col(self, size_t j):
        cdef size_t k = self._index_for_vec(1,j)
        return np.array(self._vec[k:k+self._m])

    cpdef mat(self):
        mat = np.asarray([np.asarray(self.output_col(j)) for j in get_range(1,self._n)])
        return mat

    cpdef copy(self,np.ndarray[np.float_t]  vec, size_t start_i=1,size_t start_j=1):
       cdef size_t k = self._index_for_vec(start_i,start_j)
       cdef size_t mm
       for mm in range(len(vec)):
            self._vec[k+mm]=vec[mm]

    cpdef output_as_vec(self):
        cdef C_Vec tmp = C_Vec(self._m*self._n)
        tmp.copy(self.output())
        return np.array(tmp)

    cpdef resize(self,int new_m,int new_n):

        cdef int max_m = self._m
        cdef int max_n = self._n
        if new_m < self._m or new_n < self._n:
            max_m = new_m
            max_n = new_n
        cdef np.ndarray[float] new_vec = np.zeros(new_n*new_m, self._complex)
        cdef size_t k2 = 0
        cdef size_t j, i, k
        for j in range(max_n):
            for i in range(max_m):
                k = self._index_for_vec(i+1,j+1)
                new_vec[k2] = self._vec[k]
                k2+=1
            k2+= (new_m -self._m)
        self._m = new_m
        self._n = new_n
        self._vec = new_vec

    cpdef dimensions(self):
        return 2


class Matrix_3D(Vec):
    """
    Matrix_3D is a 3 dimensional matrix.
    Under the hood this is also a vector,
    but can take 3 indicies like fortran.
    """
    def __init__(self, m,n,p, complex=False):
        super().__init__(n*m*p, complex)
        self._m = m
        self._n = n
        self._p = p

    def _index_for_vec(self, i, j, k):
        return (k-1)*(self._n*self._m)+(j-1)*self._m+i-1

    def fill(self, value, n, start_i=1,start_j=1,start_k=1):
        k = self._index_for_vec(start_i,start_j,start_k)
        super().fill(value,n,k)

    def __call__(self,i,j,k):
        return self._vec[self._index_for_vec(i,j,k)]

    def set(self, i,j, k, value):
        self._vec[self._index_for_vec(i,j,k)] = value

    def _print(self):
        print("i, j, k, value")
        for k in range(self._p):
          for j in range(self._n):
            for i in range(self._m):
                kk = self._index_for_vec(i+1,j+1,k+1)
                print(i+1,j+1,k+1,self._vec[kk])

    def get_from(self, i,j,k):
        kk = self._index_for_vec(i,j,k)
        return self._vec[kk:]

    def output_range(self, i,j,k,end):
        kk = self._index_for_vec(i,j,k)
        return self._vec[kk:kk+end]

    def copy(self, vec, start_i=1, start_j=1, start_k=1):
       k = self._index_for_vec(start_i,start_j,start_k)
       self._vec[k:k+len(vec)] = vec

    def dimensions(self):
        return 3


class Storage(object):
    """
    Fortran often used files to dump data to read back in later.
    Storage is a class to "emulate that behaviour"
    """
    def __init__(self):
        self.store = {}
        self._files = []
        self._key = {}

    def open(self, unit, file):
        if file not in self._files:
            self._files.append(file)
        self._key[unit] = file

    def close(self, unit):
        del self._key[unit]

    def write(self, unit, msg):
        if unit in self._key.keys():
           file = self._key[unit]
           if file  in self.store.keys():
               self.store[file].append(str(msg))
           else:
            self.store[file] = [str(msg)]

    def read(self, unit):
       if unit in self._key.keys():
           file = self._key[unit]
           for line in self.store[file]:
               print(line)
           print("--END--")
           print()

    def read_all(self):
        for file in self._files:
            print("---READING--- "+file)
            self.open(1, file)
            self.read(1)
            self.close(1)

    def dump(self, overwrite):
        file_path = os.path.realpath(__file__)
        write_status = "a"
        if overwrite:
            write_status="w"
        for file_name in self.store.keys():
            full_path = os.path.join(file_path, "..","..",file_name)
            print("writting to", full_path)
            new_file = open(full_path, write_status)
            for line in self.store[file_name]:
                new_file.write(line)
                new_file.write("\n")
            new_file.close()


# often want to find an index where func is true
def find_index(data,int start,int end, condition,int step =1):
    cdef int j
    for j in get_range(start,end,step):
        if condition(*data,j):
            return j
    return end


def debug_dump(file_name, val, store):
    store.open(1, file_name)
    for k in range(len(val)):
        store.write(1, val[k])
    store.close(1)


def deprecated(func):
    def wrapper(*arg):
        print("WARNING deprecated function: "+func.__name__)
        result = func(*arg)
        return result
    return wrapper
