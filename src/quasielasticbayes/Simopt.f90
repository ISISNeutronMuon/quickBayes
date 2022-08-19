      SUBROUTINE SIMOPT(X,DX,COVAR,N)
C
C Purpose
C     A SIMPLEX-based routine which optimises the values of the N 
C parameters pertaining to the components of the vector X; it also 
C estimates the related covariance matrix.
C
C Parameters
C   ARGUMENT  TYPE  I/O  DIMENSION  DESCRIPTION
C    N        I*4    I       -      No. of parameters to be optimised.
C    X        R*4    I       N      Initial guess.
C    X        R*4    O       N      The optimal answer.
C    DX       R*4    I       N      Initial step-lengths for X.
C    COVAR    R*4    O     N x N    The covariance matrix.
C
C Other Requirements
C     The user must provide a FUNCTION CCHI(X) which evalutes 
C Chi-squared (un-normalised) given the vector X.
C
C History
C     D. S. Sivia    9 Feb 1995  Initial release.
C     D. S. Sivia   21 Nov 1995  Corrected a minor, but crashable, bug.
C     D. S. Sivia   26 Apr 1996  Removed chi-squared write statements.
C
C-----------------------------------------------------------------------
C
      REAL     X(N),DX(N),COVAR(N,N)
      REAL     HESS(2500),DELTA(50),C0,C1(2,50),C2(2500)
      REAL     V(2550),EX(150),C(51)
      INTEGER  IR(51),INDX(50)
      EXTERNAL CCHI
      DATA     NMAX /50/
C
      IF (N.GT.NMAX) STOP' Sorry, too many parameters !'
      CALL VCOPY(X,V,N)
      CHIMIN=CCHI(V)
   1  CALL SIMPLEX(V,N,DX,EX,C,IR,N*1000)
      CALL VCOPY(EX,V,N)
      CHI=CCHI(V)
      IF ((1.0-CHI/CHIMIN).GT.0.0001) THEN
        CHIMIN=CHI
        GOTO 1
      ENDIF
      CALL VCOPY(V,X,N)
      CALL CHINIT(V,N,DELTA,C0,C1,C2)
      CALL HSINT1(C0,C1,DELTA,N,HESS)
      CALL HSINT2(C0,C1,C2,DELTA,N,HESS)
      CALL INVERT(HESS,COVAR,N,INDX,DETLOG)
      END
C
C***<Simplex>*****************************************************
C
      SUBROUTINE SIMPLEX(V,N,D,EX,C,IR,MX)
      REAL    V(N,*),EX(N,*),C(*),D(*)
      INTEGER IR(*)
      LOGICAL SLOW
      DATA    ALPHA,BETA,GAMA/1.0,0.5,2.0/
C
      SLOW=.FALSE.
      N3=3*N
      CAIM=0.0
      C(1)=CCHI(V(1,1))
      CALL SIMP0(V,D,C,IR,N)
      CMIN=C(IR(N+1))
      ITER=0
      NOLUCK=0
      IRSTRT=0
  10  ITER=ITER+1
      NOLUCK=NOLUCK+1
      IF (ITER.GE.MX .OR. NOLUCK.GT.100) THEN
        CALL SSORT(C,IR,N)
        CALL VCOPY(V(1,IR(N+1)),EX,N)
        IF (ITER.LT.250) THEN
          NOLUCK=0
          do ID=1,N
            D(ID)=D(ID)/10.0
          end do
          CALL VCOPY(EX,V,N)
          C(1)=CCHI(V(1,1))
          CALL SIMP0(V,D,C,IR,N)
          CMIN=C(IR(N+1))
        ELSE
          RETURN
        ENDIF
      ENDIF
      CALL VRFILL(EX,0.,N3)
      CALL XCENT(V,EX(1,3),IR(1),N)
      CALL XZERO(V(1,IR(1)),EX,EX(1,3),N,ALPHA)
      C0=CCHI(EX)
      CL=C(IR(N+1))
      CH=C(IR(1))
      CS=C(IR(2))
      IF (C0.LT.CL) THEN
        CALL EXPAND(CL,EX,EX(1,2),EX(1,3),N,C0,GAMA)
      ELSEIF (C0.GT.CS) THEN 
        CALL CONTRACT(CH,C0,C00,EX,EX(1,2),EX(1,3),V(1,IR(1)),BETA,N)
        IF (C00.LT.CH.AND.C00.LT.C0) THEN
           CALL VCOPY(EX(1,2),EX,N)
           C0=C00
        ELSE
           CALL CONT2(IR(N+1),V,C,N,IR)
           IRSTRT=IRSTRT+1
           IF (C(IR(N+1)).LT.CMIN) CMIN=C(IR(N+1))
           IF (IRSTRT.GE.5) THEN
             NOLUCK=0
             CALL SSORT(C,IR,N)
             CALL VCOPY(V(1,IR(N+1)),EX,N)
             RETURN
          ENDIF
        ENDIF
      ENDIF
C
      CALL VCOPY(EX,V(1,IR(1)),N)
      C(IR(1))=C0
      CALL SSORT(C,IR,N)
      IF (C(IR(N+1)).LT.CMIN) THEN
        DROP=(CMIN-C(IR(N+1)))/CMIN
        IF (ABS(DROP).LT.1.0E-4) SLOW=.TRUE.
        NOLUCK=0
        CMIN=C(IR(N+1))
      ENDIF
      IF (CMIN.GT.CAIM .AND. .NOT. SLOW) GOTO 10
      CALL VCOPY(V(1,IR(N+1)),EX,N)
      END
C     ----------------------------
      SUBROUTINE SIMP0(V,D,C,IR,N)
      REAL    V(N,*),C(*),D(*)
      INTEGER IR(*)
      do I=2,N+1
        CALL VCOPY(V,V(1,I),N)
        V(I-1,I)=V(I-1,I)+D(I-1)
        C(I)= CCHI(V(1,I))
      end do
      CALL SSORT(C,IR,N)
      END
C     ------------------------
      SUBROUTINE SSORT(C,IR,N)
      REAL    C(*)
      INTEGER IR(*)
      IR(1)=1
      DO 20 J=2,N+1
        DO 10 I=1,J-1
          IF (C(J).LT.C(IR(I))) GOTO 10
          do II=1,J-I
            IR(J+1-II)=IR(J-II)
          end do
          IR(I)=J
          GOTO 20
  10    CONTINUE
        IR(J)=J
  20  CONTINUE
      END
C     ---------------------------
      SUBROUTINE XCENT(V,XC,IH,N)
      REAL V(N,*),XC(*)
      XNORM=1.0/FLOAT(N)
      DO 20 J=1,N+1
        IF (J.EQ.IH) GOTO 20
        do I=1,N
          XC(I)=XC(I)+V(I,J)
        end do
  20  CONTINUE
      do I=1,N
        XC(I)=XC(I)*XNORM
      end do
      END
C     ------------------------------
      SUBROUTINE XZERO(XH,X0,XC,N,A)
      REAL XH(*),X0(*),XC(*)
      do I=1,N
        X0(I)=A*(XC(I)-XH(I))+XC(I)
      end do
      END
C     --------------------------------------
      SUBROUTINE EXPAND(CL,X0,X00,XC,N,C0,G)
      REAL X0(*),X00(*),XC(*)
      do I=1,N
        X00(I)=G*(X0(I)-XC(I))+XC(I)
      end do
      C00=CCHI(X00)
      IF (C00.LT.CL) THEN
        CALL VCOPY(X00,X0,N)
        C0=C00
      ENDIF
      END
C     --------------------------------------------------
      SUBROUTINE CONTRACT(CH,C0,C00,X0,X00,XC,XH,B,N)
      REAL X0(*),X00(*),XC(*),XH(*)
      IF (C0.LT.CH) THEN
        do I=1,N
          X00(I)=B*(X0(I)-XC(I))+XC(I)
        end do
      ELSE
        do I=1,N
          X00(I)=B*(XH(I)-XC(I))+XC(I)
        end do
      ENDIF
      C00=CCHI(X00)
      END
C     ----------------------------
      SUBROUTINE CONT2(K,V,C,N,IR)
      REAL     V(N,*),C(*)
      INTEGER  IR(*)
      DO 20 J=1,N+1
        IF (J.EQ.K) GOTO 20
        do I=1,N
          V(I,J)=0.5*(V(I,J)+V(I,K))
          C(J)=CCHI(V(1,J))
        end do
  20  CONTINUE
      CALL SSORT(C,IR,N)
      END
C
C***<error analysis>****************************************************
C
      SUBROUTINE CHINIT(V,N,DELTA,C0,C1,C2)
      REAL     V(*),DELTA(*),C0,C1(2,*),C2(N,*)
      EXTERNAL CCHI
      CALL VRFILL(C2,0.0,N*N)
      C0=CCHI(V)
      CDELTA=C0/2000.0
      IF (CDELTA.LT.0.5) CDELTA=0.5
      do I=1,N
        CALL DLINIT(V,I,C1(1,I),C1(2,I),C0,CDELTA,DELTA(I))
      end do
      do J=1,N
        V(J)=V(J)+DELTA(J)
        do I=J+1,N
          V(I)=V(I)+DELTA(I)
          C2(I,J)=CCHI(V)
          C2(J,I)=C2(I,J)
          V(I)=V(I)-DELTA(I)
        end do
        V(J)=V(J)-DELTA(J)
      end do
      END
C     ----------------------------------------------------
      SUBROUTINE DLINIT(X,J,CMINUS,CPLUS,C0,CDELTA,DELTA)
      REAL X(*)
      DATA FRAC,SMALL,XLARGE /1.0E-6,1.0E-20,1.0E20/
      X0=X(J)
      D=FRAC*ABS(X0)+SMALL
      do I=1,100
        D=D+D
        X(J)=X0-D
        CMINUS=CCHI(X)
        X(J)=X0+D
        CPLUS=CCHI(X)
        IF (X(J).GT.XLARGE) GOTO 1
        IF (ABS(CMINUS-C0).GT.CDELTA) GOTO 1
        IF (ABS(CPLUS-C0).GT.CDELTA) GOTO 1
      end do
   1  DELTA=D
      X(J)=X0
      END
C     -------------------------------
      SUBROUTINE HSINT1(C0,C1,D,N,HS)
      REAL C0,C1(2,*),D(*),HS(N,*)
      do I=1,N
        HS(I,I)=(C1(2,I)-C0-C0+C1(1,I))/(D(I)*D(I))
      end do
      END
C     ----------------------------------
      SUBROUTINE HSINT2(C0,C1,C2,D,N,HS)
      REAL C0,C1(2,*),C2(N,*),D(*),HS(N,*)
      do J=1,N
        do I=J+1,N
          HS(I,J)=(C2(I,J)-C1(2,I)-C1(2,J)+C0)/(D(I)*D(J))
          HS(J,I)=HS(I,J)
        end do
      end do
      END
