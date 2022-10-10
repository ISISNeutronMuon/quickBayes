from fortran_python import get_range, complex_zeros
import numpy as np
import scipy.fft as sc

import cython
from libc.math cimport sin, pi
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# It's necessary to call "import_array" if you use any part of the
# numpy PyArray_* API. From Cython 3, accessing attributes like
# ".shape" on a typed Numpy array use this API. Therefore we recommend
# always calling "import_array" whenever you "cimport numpy"
np.import_array()


def CMPLX(A, B):
    return complex(A, B)


def Re(x):
    y = []
    for k in x:
        y.append(k.real)
    return np.asarray(y)


def Im(x):
    y = []
    for k in x:
        y.append(k.imag)
    return np.asarray(y)


def flatten_cmplx(np.ndarray[np.complex_t] x):
    cdef int N = len(x)
    cdef int v

    cdef np.ndarray[np.float_t] y = np.zeros(N * 2)
    for v in range(N):
        y[2 * v] = x[v].real
        y[2 * v + 1] = x[v].imag
    return y


def flatten_real(np.ndarray[np.float_t] x):
    cdef int N = len(x)
    cdef int v

    cdef np.ndarray[np.float_t] y = np.zeros(N * 2)
    for v in range(N):
        y[2 * v] = x[v]
        y[2 * v + 1] = 0.0
    return y


def flatten(x):
    if any(np.iscomplex(x)):
        return flatten_cmplx(x)
    else:
        return flatten_real(x)


def compress(np.ndarray[np.float_t] x):
    y = []
    cdef int end = len(x)
    cdef int k
    if end % 2 != 0:
        end -= 1
    for k in range(0, end, 2):
        y.append(complex(x[k], x[k + 1]))
    return np.asarray(y)


def FOUR2_NEG_IFORM(np.ndarray[np.float_t] DATA, int N, int NDIM, int ISIGN, int IFORM):
    # NDIM is always 1
    cdef int N_tot = 1
    N_tot = N_tot * N
    N_tot = (N_tot / N) * (N / 2 + 1)
    cdef int NREM = 1
    cdef int JDIM = 1
    cdef int IDIM = 1
    cdef int NCURR = N
    NCURR = int(N / 2)
    cdef int nn = N
    DATA = FOXRL(DATA, nn, NREM, ISIGN, IFORM)
    N_tot = N_tot / (N / 2 + 1) * N
    cdef int NPREV = N_tot / (N * NREM)
    DATA = flatten(BOTRV(compress(DATA), NPREV, NCURR, NREM))
    cdef np.ndarray[np.complex_t] result = COOL2(compress(DATA), NPREV, NCURR, NREM, ISIGN)
    return result


def FOXRL(np.ndarray[np.float_t] DATA, int N, int NREM, int ISIGN, int IFORM):
    cdef float TWOPI = 2 * pi * float(ISIGN)
    cdef int IP0 = 2
    cdef int IP1 = IP0 * int(float(N) / 2.)
    cdef int IP2 = IP1 * NREM
    cdef float TEMPR, THETA, SINTH, ZSTPR, ZSTPI, ZR, ZI, DIFR, DIFI, TEMPI
    cdef int I1MIN, I1MAX, I1, I2, I2CNJ
    cdef int J1 = IP1 + 1
    DATA[1] = DATA[J1 - 1]

    for I2 in get_range(1, IP2, IP1):
        TEMPR = DATA[I2 - 1]
        DATA[I2 - 1] = DATA[I2 - 1] + DATA[I2]
        DATA[I2] = TEMPR - DATA[I2]
    THETA = TWOPI / float(N)
    SINTH = sin(THETA / 2.)
    ZSTPR = -2. * SINTH * SINTH
    ZSTPI = sin(THETA)
    ZR = (1. - ZSTPI) / 2.
    ZI = (1. + ZSTPR) / 2.
    ZR = 1. - ZR
    ZI = -ZI
    I1MIN = IP0 + 1
    I1MAX = IP0 * int(N / 4) + 1
    for I1 in get_range(I1MIN, I1MAX, IP0):

        for I2 in get_range(I1, IP2, IP1):

            I2CNJ = int(IP0 * (N / 2 + 1) - 2 * I1 + I2)
            if I2 - I2CNJ >= 0:
                if ISIGN * (2 * IFORM + 1) < 0:
                    DATA[I2] = -DATA[I2]
            else:
                DIFR = DATA[I2 - 1] - DATA[I2CNJ - 1]
                DIFI = DATA[I2] + DATA[I2CNJ]
                TEMPR = DIFR * ZR - DIFI * ZI
                TEMPI = DIFR * ZI + DIFI * ZR
                DATA[I2 - 1] = DATA[I2 - 1] - TEMPR
                DATA[I2] = DATA[I2] - TEMPI
                DATA[I2CNJ - 1] = DATA[I2CNJ - 1] + TEMPR
                DATA[I2CNJ] = DATA[I2CNJ] - TEMPI
                DATA[I2CNJ - 1] = 2 * DATA[I2CNJ - 1]
                DATA[I2CNJ] = 2 * DATA[I2CNJ]
            DATA[I2 - 1] = 2 * DATA[I2 - 1]
            DATA[I2] = 2 * DATA[I2]

        TEMPR = ZR - .5
        ZR = ZSTPR * TEMPR - ZSTPI * ZI + ZR
        ZI = ZSTPR * ZI + ZSTPI * TEMPR + ZI

    return DATA


def BOTRV(np.ndarray[np.complex_t] DATA, int NPREV, int N, int NREM):
    cdef int IP0 = 1
    cdef int IP1 = IP0 * NPREV
    cdef int IP4 = IP1 * N
    cdef int IP5 = IP4 * NREM
    cdef int I4REV = 1
    cdef int I4MAX = IP4
    cdef int I1MAX = 1
    cdef int IP2 = 1
    cdef np.complex_t TEMP

    cdef int I4, I5REV, I5
    for I4 in get_range(1, I4MAX, IP1):
        if I4 - I4REV < 0:
            I1MAX = I4 + IP1 - IP0
            for I1 in get_range(I4, I1MAX, IP0):
                for I5 in get_range(I1, IP5, IP4):
                    I5REV = I4REV + I5 - I4
                    TEMP = DATA[I5 - 1]
                    DATA[I5 - 1] = DATA[I5REV - 1]
                    DATA[I5REV - 1] = TEMP
        IP2 = int(IP4 / 2)
        if I4REV - IP2 > 0:
            while IP2 - IP1 > -1:
                if I4REV - IP2 <= 0:
                    break
                I4REV = I4REV - IP2
                IP2 = int(IP2 / 2)
        I4REV = I4REV + IP2
    return DATA



def COOL2(np.ndarray[double complex] DATA, int NPREV, int N, int NREM, int ISIGN):

    # can't do the complex_t's in this function as it causes arithmatic errors
    cdef float TWOPI = 2 * pi * float(ISIGN)
    cdef int IP0 = 1
    cdef int IP1 = IP0 * NPREV
    cdef int IP4 = IP1 * N
    cdef int IP5 = IP4 * NREM
    cdef int IP2 = IP1
    cdef int NPART = N
    cdef int IP3 = 0
    flag = True
    cdef int I1, I1Max, I3A, I3B, I3C, I3D, I5, I2
    cdef double complex  TMP, T0, T1, T2, TEMP, WSTP, W, W2, W3
    cdef float THETA, SINTH
    while NPART - 2 > 0:
        NPART = int(NPART / 4)

    flag = True
    if NPART - 2 == 0:
        IP3 = IP2 * 2
        I1MAX = IP1
        for I1 in get_range(1, I1MAX, IP0):
            for I5 in get_range(I1, IP5, IP3):
                I3A = I5
                I3B = I3A + IP2
                TMP = DATA[I3B - 1]
                DATA[I3B - 1] = DATA[I3A - 1] - TMP
                DATA[I3A - 1] = DATA[I3A - 1] + TMP
                flag = True
        IP2 = IP3
        if IP2 - IP4 >= 0:
            return DATA
    while IP3 - IP4 < 0:

        IP3 = IP2 * 4
        THETA = TWOPI / float(IP3 / IP1)
        SINTH = sin(THETA / 2.)
        WSTP = np.cdouble(complex(-2. * SINTH * SINTH, sin(THETA)))
        W = 1. + 0 * 1j
        for I2 in get_range(1, IP2, IP1):
            if I2 - 1 > 0:
                W2 = pow(W, 2) + 0 * 1j
                W3 = pow(W, 3) + 0 * 1j
            else:
                W2 = 1. + 0 * 1j
                W3 = 1. + 0 * 1j
            I1MAX = I2 + IP1 - IP0
            for I1 in get_range(I2, I1MAX, IP0):
                for I5 in get_range(I1, IP5, IP3):
                    I3A = I5
                    I3B = int(I3A + IP2)
                    I3C = int(I3B + IP2)
                    I3D = int(I3C + IP2)
                    if I2 - 1 > 0:
                        DATA[I3B - 1] = W2 * DATA[I3B - 1]
                        DATA[I3C - 1] = W * DATA[I3C - 1]
                        DATA[I3D - 1] = W3 * DATA[I3D - 1]
                    T0 = DATA[I3A - 1] + DATA[I3B - 1]
                    T1 = DATA[I3A - 1] - DATA[I3B - 1]
                    T2 = DATA[I3C - 1] + DATA[I3D - 1]
                    T3 = DATA[I3C - 1] - DATA[I3D - 1]
                    DATA[I3A - 1] = T0 + T2
                    DATA[I3C - 1] = T0 - T2
                    TEMP = 1j * T3
                    if ISIGN < 1:
                        TEMP = -TEMP
                    DATA[I3B - 1] = T1 + TEMP
                    DATA[I3D - 1] = T1 - TEMP
            W = W * WSTP + W
        IP2 = IP3
    return DATA


def call_scipy(data):
    result = sc.fft(data)
    return result


def do_transform(np.ndarray[np.float_t] tmp, int N, int NDIM, int ISIGN, int IFORM):
    cdef np.ndarray[np.complex_t] result
    if IFORM == -1:
        """
        N < len(DATA) since we dont actually use all of the data
        This seems to take data of the form exp(-t*gamma)*cos(omega*t) and
        pick out the omega freq -> the oscialltion freq is related to the
        mean of the gaussian. The envelope is related to the width of the
        gaussianss
        """
        result = FOUR2_NEG_IFORM(tmp, N, NDIM, ISIGN, IFORM)
    elif ISIGN == 1 and IFORM == 1:  # assume real data
        result = np.conj(call_scipy(tmp))
    elif ISIGN == -1 and IFORM == 1:
        result = call_scipy(tmp)
    elif ISIGN == 1 and IFORM == 0:
        result = np.conj(call_scipy(tmp[0:N]))
        result = result[0:int(N / 2) + 1]
    elif ISIGN == -1 and IFORM == 0:
        result = call_scipy(tmp[0:N])
        result = result[0:int(N / 2) + 1]
    return result


def FOUR2(DATA, N, NDIM, ISIGN, IFORM):
    tmp = DATA.output()
    if any(np.iscomplex(tmp)):
        tmp = flatten(tmp)
    return do_transform(tmp, N, NDIM, ISIGN, IFORM)


def FOUR2_raw(np.ndarray[np.float_t] DATA, int N, int NDIM, int ISIGN, int IFORM):
    # we know all of the calls have IFORM = -1
    return FOUR2_NEG_IFORM(DATA, N, NDIM, ISIGN, IFORM)


def FOUR2_IFT(data, int N, int NDIM, int ISIGN):
    if any(np.iscomplex(data)):
        data = flatten(data)

    # This seems to take data of the form exp(-t*gamma)*cos(omega*t) and pick
    # out the omega freq -> the oscialltion freq is related to the mean of the
    # gaussian. The envelope is related to the width of the gaussianss
    result = FOUR2_NEG_IFORM(data, N, NDIM, ISIGN, -1)
    return result


def FOUR3(DATA, N, NDIM, ISIGN, IFORM):
    return np.asarray(sc.ifft(DATA.output()[0:N]))


def four4(DATA):
    return np.asarray(sc.ifft(DATA))
