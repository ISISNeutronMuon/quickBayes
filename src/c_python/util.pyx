from quasielasticbayes.fortran_python import Vec, get_range, round_sig
import math

"""
These could probably be made faster
These (supposedly) do Lower matrix decomposition.
Cannot get a match with scipy
"""


def LUDCMP(A, N, NP):
    """
    To match the fortran with the test data (quasi-lines =2)
    we need to round the values to a lower level of precision.
    Otherwise the errors become huge.
    """
    INDX = Vec(N)
    NMAX = 100
    TINY = 1.0E-20
    VV = Vec(N)
    DUM = 0
    IMAX = 0
    D = 1.0

    for I in get_range(1, N):
        AAMAX = 0.0
        for J in get_range(1, N):
            if round_sig(abs(A(I, J)), 8) > AAMAX:
                AAMAX = round_sig(abs(A(I, J)), 8)

        if AAMAX == 0.0:
            print(' Singular matrix!')
            return
        VV.set(I, round_sig(1.0 / AAMAX, 8))

    for J in get_range(1, N):
        if J > 1:
            for I in get_range(1, J - 1):
                SUM = round_sig(A(I, J), 8)
                if I > 1:
                    for K in get_range(1, I - 1):
                        SUM = round_sig(SUM -
                                        round_sig(A(I, K), 8) *
                                        round_sig(A(K, J), 8))
                    A.set(I, J, SUM)

        AAMAX = 0.0
        for I in get_range(J, N):
            SUM = round_sig(A(I, J), 8)
            if J > 1:
                for K in get_range(1, J - 1):
                    SUM = round_sig(SUM - round_sig(A(I, K), 8)
                                    * round_sig(A(K, J), 8), 8)
                A.set(I, J, SUM)

            DUM = round_sig(round_sig(VV(I), 8) * abs(round_sig(SUM, 8)), 8)
            if round_sig(DUM, 8) >= round_sig(AAMAX, 8):
                IMAX = I
                AAMAX = round_sig(DUM, 8)

        if J != IMAX:
            for K in get_range(1, N):
                DUM = round_sig(A(IMAX, K), 8)
                A.set(IMAX, K, round_sig(A(J, K), 8))
                A.set(J, K, DUM)

            D = -D
            VV.set(IMAX, round_sig(VV(J), 8))

        INDX.set(J, int(IMAX))
        if J != N:
            if A(J, J) == 0.0:
                A.set(J, J, TINY)
            DUM = 1.0 / A(J, J)
            for I in get_range(J + 1, N):
                A.set(I, J, round_sig(A(I, J) * DUM, 8))

    if A(N, N) == 0.0:
        A.set(N, N, TINY)
    return INDX, D, A


def LUBKSB(A, N, NP, INDX, B):
    # B is a numpy array
    II = 0
    for I in get_range(1, N):
        LL = int(INDX(I))
        SUM = B[LL - 1]
        B[LL - 1] = B[I - 1]
        if II != 0:
            for J in get_range(II, I - 1):
                SUM = SUM - A(I, J) * B[J - 1]

        elif SUM != 0.0:
            II = I

        B[I - 1] = SUM
    for I in get_range(N, 1, -1):
        SUM = B[I - 1]
        if I < N:
            for J in get_range(I + 1, N):
                SUM = SUM - A(I, J) * B[J - 1]

        B[I - 1] = SUM / A(I, I)
    return B
