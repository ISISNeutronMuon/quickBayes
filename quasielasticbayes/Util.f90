      SUBROUTINE open_f(nit,fname)
      character*(*) fname
      logical found
      integer nit
      INQUIRE(FILE=fname,EXIST=found)
      if(found)then
       OPEN(UNIT=nit,FILE=fname,STATUS='OLD',FORM='FORMATTED')
       CLOSE(UNIT=nit,STATUS='DELETE')
      endif
      OPEN(UNIT=nit,FILE=fname,STATUS='NEW',FORM='FORMATTED')
      RETURN
      END

C***<initialise arrays>*************************************************
C
      SUBROUTINE VCOPY(X,Y,N)
      REAL X(*),Y(*)
      do I=1,N
        Y(I)=X(I)
      end do
      END
C     --------------------------
      SUBROUTINE VIFILL(IN,K,N)
      INTEGER IN(*)
      do I=1,N
        IN(I)=K
      end do
      END
C     ------------------------
      SUBROUTINE VRFILL(X,A,N)
      REAL X(*)
      do I=1,N
        X(I)=A
      end do
      END
C     ------------------------
      SUBROUTINE VLFILL(LG,L,N)
      LOGICAL LG(*),L
      do I=1,N
        LG(I)=L
      end do
      END
C----------------------------------------------------------------------C
C  Two subroutines, for Cholesky decomposition of a matrix, copied out C
C  of Numerical Recipes.                                               C
C----------------------------------------------------------------------C
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION  A(NP,NP),INDX(N),VV(NMAX)
      D=1.0
      do I=1,N
        AAMAX=0.0
        do J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
        end do
        IF (AAMAX.EQ.0.0) STOP ' Singular matrix!'
        VV(I)=1.0/AAMAX
      end do
      do J=1,N
        IF (J.GT.1) THEN
          do I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1) THEN
              do K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
              end do
              A(I,J)=SUM
            ENDIF
          end do
        ENDIF
        AAMAX=0.0
        do I=J,N
          SUM=A(I,J)
          IF (J.GT.1) THEN
            do K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
            end do
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
        end do
        IF (J.NE.IMAX) THEN
          do K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
          end do
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF (J.NE.N) THEN
          IF (A(J,J).EQ.0.0) A(J,J)=TINY
          DUM=1.0/A(J,J)
            do I=J+1,N
              A(I,J)=A(I,J)*DUM
            end do
        ENDIF
      end do
      IF (A(N,N).EQ.0.0) A(N,N)=TINY
      END
C     --------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      do I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0) THEN
          do J=II,I-1
            SUM=SUM-A(I,J)*B(J)
          end do
        ELSEIF (SUM.NE.0.0) THEN
          II=I
        ENDIF
        B(I)=SUM
      end do
      do I=N,1,-1
        SUM=B(I)
        IF (I.LT.N) THEN
          do J=I+1,N
            SUM=SUM-A(I,J)*B(J)
          end do
        ENDIF
        B(I)=SUM/A(I,I)
      end do
      END
