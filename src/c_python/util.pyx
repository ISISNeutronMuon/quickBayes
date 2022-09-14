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


def test(a, n, np, indx, d):
    d = 1
    # looking for the largest a in each row and store it in vv as inverse
    # We need a new list same size as indx, for this we use .copy()
    vv = indx.copy()
    for i in range(0, n):
        big = 0.0
        for j in range(0, n):
            temp = math.fabs(a[i][j])
            if (temp > big):
                big = temp
        vv[i] = 1.0 / big
    #
    # run Crout's algorithm
    for j in range(0, n):
        # top half & bottom part are combined
        # but the upper limit l for k sum is different
        big = 0.0
        for i in range(0, n):
            if (i < j):
                l = i
            else:
                l = j
            sum = a[i][j]
            for k in range(0, l):
                sum -= a[i][k] * a[k][j]
            a[i][j] = sum
            # for bottom half, we keep track which row is larger
            if (i >= j):
                dum = vv[i] * math.fabs(sum)
                if (dum >= big):
                    big = dum
                    imax = i
        # pivoting part, swap row j with row imax, a[j] is a whole row
        if (j != imax):
            dum = a[imax]
            a[imax] = a[j]
            a[j] = dum
            d = - d
            vv[imax] = vv[j]
        # divide by the beta diagonal value
        indx[j] = int(imax)
        dum = 1.0 / a[j][j]
        for i in range(j + 1, n):
            a[i][j] *= dum
    return a, indx, d


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
