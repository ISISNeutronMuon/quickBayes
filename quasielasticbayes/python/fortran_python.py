# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

import numpy as np
import math
from scipy import linalg
import os
"""
This file contains methods and classes to make converting Fortran to Python easier.
This is achived by this class making adjustments so that the syntax of the original
code is changed by as little as possible
"""

def get_range(start, end, dx=1):
    return range(start, int(end+1*np.sign(dx)),int(dx)) 

def NINT(value):
    return int(round(value))

def CMPLX(A,B):
    return np.cdouble(complex(A,B))

def complex_zeros(n):
    return np.asarray([complex(0,0) for _ in range(n)])

def zeros(n, complex=False):
    if complex:
        return complex_zeros(n)
    else:
        return np.zeros(n)

class vec(object):
    def __init__(self, n,complex=False):
        self._vec = zeros(n, complex)
        
    def __call__(self,j):
        return self._vec[j-1]
        
    def _print(self):
        print(self._vec)

    def testing(self,j):
        print("testing",j)
        
    def set(self, j, value):
        self._vec[j-1] = value

    def output(self):
        return self._vec

    def output_range(self,start=1,end=2):
        return self._vec[start-1:end]

    def fill(self, value, n,start=1):
        self._vec[start-1:n+start-2] = value
        
    def get_from(self, start=1):
        return self._vec[start-1:]
        
    def copy(self, vector, start_i=1):
       self._vec[start_i-1:start_i-1+len(vector)] = vector

    #def update(self,vector, start):
    #    self._vec[start-1:] = vector

    def dimensions(self):
        return 1

class BoolVec(vec):
    def __init__(self, n,complex=False):
        self._vec = [True for _ in range(n)]

class matrix_2(vec):
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
   
                          
    def dimensions(self):
        return 2
        
        
class matrix_3(vec):
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
    
    def __call__(self,i,j):
        return self._vec[self._index_for_vec(i,j)]
        
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
    
    def copy(self, vec, start_i=1, start_j=1, start_k=1):
       k = self._index_for_vec(start_i,start_j,start_k)
       self._vec[k:k+len(vec)] = vec  
                          
    def dimensions(self):
        return 3    
        
"""
Fortran often used files to dump data to read back in later.
This is a class to "emulate that behaviour"
"""
class storage(object):
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

    def dump(self):
        file_path = os.path.realpath(__file__)
        for file_name in self.store.keys():
            full_path = os.path.join(file_path, "..","..",file_name)
            print("writting to", full_path)
            new_file = open(full_path, "w")
            for line in self.store[file_name]:
                new_file.write(line)
                new_file.write("\n")
            new_file.close()

# often want to find an index where func is true
def find_index(data, start, end, condition, step =1):
    for j in get_range(start,end,step):
        if condition(*data,j):
            return j
    return end
