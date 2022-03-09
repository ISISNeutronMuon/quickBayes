# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
from  quasielasticbayes.python.fortran_python import *


"""
These could probably be made faster
These (supposedly) do Lower matrix decomposition. 
Cannot get a match with scipy
"""

def LUDCMP(A,N,NP):
      NMAX=100
      TINY=1.0E-20
      VV = vec(NMAX)
      DUM=0
      IMAX=0
      INDX = vec(N)
      #DIMENSION  A(NP,NP),INDX(N),VV(NMAX)
      D=1.0
      for I in get_range(1,N):
        AAMAX=0.0
        for J in get_range(1,N):
          if abs(A(I,J)) > AAMAX:
              AAMAX=abs(A(I,J))
        
        if AAMAX == 0.0:
            print(' Singular matrix!')
            return
        VV.set(I,1.0/AAMAX)
      for J in get_range(1,N):
        if J > 1:
          for I in get_range(1,J-1):
            SUM=A(I,J)
            if I>1:
              for K in get_range(1,I-1):
               SUM=SUM-A(I,K)*A(K,J)
              
              A.set(I,J,SUM)

        AAMAX=0.0
        for I in get_range(J,N):
          SUM=A(I,J)
          if J>1:
            for K in get_range(1,J-1):
              SUM=SUM-A(I,K)*A(K,J)
            
            A.set(I,J,SUM)
          
          DUM=VV(I)*abs(SUM)
          if DUM>=AAMAX:
            IMAX=I
            AAMAX=DUM

        if J!=IMAX:
          for K in get_range(1,N):
            DUM=A(IMAX,K)
            A.set(IMAX,K, A(J,K))
            A.set(J,K,DUM)
          
          D=-D
          VV.set(IMAX,VV(J))
        
        INDX.set(J,IMAX)
        if J!=N:
          if A(J,J)==0.0:
              A.set(J,J,TINY)
          DUM=1.0/A(J,J)
          for I in get_range(J+1,N):
              A.set(I,J,A(I,J)*DUM)

      if A(N,N)==0.0:
          A.set(N,N,TINY)
      return INDX, D

def LUBKSB(A,N,NP,INDX,B):
      #DIMENSION A(NP,NP),INDX(N),B(N)
      # B is a numpy array
      II=0
      for I in get_range(1,N):
        LL=int(INDX(I))
        SUM=B(LL)
        B.set(LL, B(I))
        if II!=0:
          for J in get_range(II,I-1):
            SUM=SUM-A(I,J)*B(J)
          
        elif SUM!=0.0:
          II=I
        
        B.set(I,SUM)
      for I in get_range(N,1,-1):
        SUM=B(I)
        if I<N:
          for J in get_range(I+1,N):
            SUM=SUM-A(I,J)*B(J)

        B.set(I,SUM/A(I,I))
