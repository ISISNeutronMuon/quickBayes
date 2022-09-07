# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
from fortran_python import round_sig

from numpy import zeros
cimport numpy as np
# It's necessary to call "import_array" if you use any part of the
# numpy PyArray_* API. From Cython 3, accessing attributes like
# ".shape" on a typed Numpy array use this API. Therefore we recommend
# always calling "import_array" whenever you "cimport numpy"
np.import_array()

# shift bin values onto new grid
def bin_shift_vecs(np.ndarray[np.float_t] y_grid,int N, np.ndarray[np.float_t] index_vec, np.ndarray[np.float_t] XPDAT): # is this the slow down?
      cdef np.ndarray[np.float_t] y_shifted = zeros(N)
      cdef int I, J
      cdef float fractional_x_shift

      for I in range(N):
        J=int(index_vec[I]) # get the index that says where the shift value is
        fractional_x_shift = XPDAT[I]
        # get fractions of original bins in the new shifted bin add sum (e.g fractional_x_shift = 0.2)
        y_shifted[I] = y_grid[J-1]*(1-fractional_x_shift) + fractional_x_shift*y_grid[J]

      return y_shifted

def HESS0_calc(np.ndarray[np.float_t] RESID, np.ndarray[np.float_t] DDDPAR, int N):
      cdef float SM = 0.0
      cdef int magic = 6
      cdef int kk
      for kk in range(N):
          SM = round_sig(SM,magic) + round_sig(round_sig(RESID[kk],magic)*round_sig(DDDPAR[kk],magic),magic)
      return SM
