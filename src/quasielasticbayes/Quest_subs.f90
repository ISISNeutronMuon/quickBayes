      SUBROUTINE PLNORM(ZLMAX,ZL,Z,N)
      REAL Z(*),ZL(*)
      do I=1,N
        ZL(I)=ZL(I)-ZLMAX
        Z(I)=EXP(ZL(I))
      end do
      END
C     ------------------------------
      SUBROUTINE PRMARG(Z,NX,NY,X,Y)
      REAL Z(NX,NY),X(*),Y(*)
      CALL VRFILL(X,0.0,NX)
      CALL VRFILL(Y,0.0,NY)
      YMAX=-1.0E10
      do J=1,NY
        SUM=0.0
        do I=1,NX
          ZIJ=Z(I,J)
          SUM=SUM+ZIJ
          X(I)=X(I)+ZIJ
        end do
        Y(J)=SUM
        YMAX=MAX(YMAX,SUM)
      end do
      XMAX=-1.0E10
      do I=1,NX
        XMAX=MAX(XMAX,X(I))
      end do
      XSCL=1.0/XMAX
      do I=1,NX
        X(I)=X(I)*XSCL
      end do
      YSCL=1.0/YMAX
      do J=1,NY
        Y(J)=Y(J)*YSCL
      end do
      END
C     ------------------------------
      SUBROUTINE init_paras
      INCLUDE 'res_par.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      COMMON /WRKCOM/ WORK(m_d2,2)
      do n=1,m_d
       XJ(n)=0.0
       XDAT(n)=0.0
       DAT(n)=0.0
       SIG(n)=0.0
       IPDAT(n)=0
       XPDAT(n)=0.0
       FIT(n)=0.0
       RESID(n)=0.0
       do m=1,m_p
        DDDPAR(n,m)=0.0
       end do
      end do
      do n=1,m_d1
       TWOPIK(n)=0.0
       do m=1,3
        EXPF(n,m)=0.0
       end do
      end do
      do n=1,m_d2
       FRES(n)=0.0
       FWRK(n)=0.0
       do m=1,2
        FR2PIK(n,m)=0.0
        WORK(n,m)=0.0
       end do
      end do
      do n=1,m_p
       FITP(n)=0.0
       do m=1,2
        SCLVEC(n,m)=0.0
       end do
      end do
      do n=1,m_sp
       QW1(n)=0.0
       SIGQW1(n)=0.0
      end do
      RETURN
      END
C
C**<some initialisation routines>***************************************
C
      SUBROUTINE BTINIT(X,N,XMIN,XMAX,NBMIN)
      REAL    X(*)
      INTEGER NBMIN(*)
      NBMIN(1)=N
      NBMIN(2)=NINT(0.1*FLOAT(N))
      DX=(XMAX-XMIN)/FLOAT(N-1)
      X(1)=XMIN
      do I=2,N
        X(I)=X(I-1)+DX
      end do
      END
C     ----------------------------------------
      SUBROUTINE SGINIT(Y,NMAX,N,Y0,YSIG,WSCL)
      REAL Y(NMAX,*)
      YMAX=Y0+MIN(20.0*YSIG,Y0)
      YMIN=YMAX/FLOAT(N)
      DY=(YMAX-YMIN)/FLOAT(N-1)
      Y(1,1)=YMIN-DY/2.0
      Y(1,2)=2.0*Y(1,1)*WSCL
      do I=2,N
        Y(I,1)=Y(I-1,1)+DY
        Y(I,2)=2.0*Y(I,1)*WSCL
      end do
      END
C
C***<calculate data & chi-squared>**************************************
C
      FUNCTION CCHI(V)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      COMMON /WRKCOM/ WORK(m_d2,2)
      REAL V(*)
      CHI=0.0
      B1=BSCL*V(1)
      B2=BSCL*V(2)
      A0=ASCL*V(3)
      DELTAX=V(4)
      NFT2=NFFT/2+1
      CALL CXSHFT(FRES,DELTAX,TWOPIK,FR2PIK,FR2PIK(1,2),NFT2)
      CALL VRFILL(WORK,A0,NFT2)
      XNSCL=-LOG(1.0E-7)/(TWOPIK(2)-TWOPIK(1))
      do J=1,NFEW
        AJ=ASCL*V(3+J+J)
        SIGJ=WSCL*V(4+J+J)/GSCL
        CALL VRFILL(EXPF(1,J),0.0,NFT2)
        NXFT2=1+NINT(XNSCL/(ABS(SIGJ)+1.0E-10))
        IF (NXFT2.GT.NFT2) NXFT2=NFT2
        do I=1,NXFT2
          EXPIJ=EXP(-TWOPIK(I)*SIGJ)
          EXPF(I,J)=EXPIJ
          WORK(I,1)=WORK(I,1)+AJ*EXPIJ
        end do
      end do
      CALL VMLTRC(WORK,FR2PIK,NFT2,FWRK)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,FIT)
      X1=XDAT(1)
      if(o_bgd.eq.2) BNRM=(B2-B1)/(XDAT(NDAT)-X1)
C     avoid conflict BNORM with ModPars
      do I=1,NDAT
        FIT(I)=FIT(I)+B1
         if(o_bgd.eq.2) FIT(I)=FIT(I)+BNRM*(XDAT(I)-X1)
        DIF=FIT(I)-DAT(I)
        RESID(I)=DIF*SIG(I)
        CHI=CHI+DIF*RESID(I)
      end do
      if(o_w1.eq.1)then
       if(NFEW.GE.1) then
        RESW1D=(WSCL*V(6)-QW1(ISPEC))/SIGQW1(ISPEC)
        CHI=CHI+2.0*RESW1D**2
       endif
      endif
      CCHI=CHI/(2.0*FLOAT(NDAT))
      END
C
C**<search for one more & refine amplitudes>****************************
C
      SUBROUTINE SEARCH(GRAD,HESS,DPAR,NFEW,INDX,COVAR,FITP)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      REAL            GRAD(*),HESS(*),DPAR(*),COVAR(*),FITP(*)
      INTEGER         INDX(*)
      J=4+2*NFEW
      DXLOG=0.85
      NSRCH=NINT(LOG(5.0*GSCL/WSCL)/LOG(DXLOG))
      CMIN=1.0E20
      FITP(J-1)=0.1
      FITP(J)=1.0
      do I=1,NSRCH
        CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
        CNORM=CBCHI(FITP)
        IF (CNORM.LT.CMIN) THEN
          CMIN=CNORM
          SIGJ=FITP(J)
        ENDIF
        FITP(J)=FITP(J)*DXLOG
      end do
      FITP(J)=SIGJ
      CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
      END
C
C**<refinement>*********************************************************
C
      SUBROUTINE REFINE(GRAD,HESS,NP,DETLOG,INDX,COVAR,STEPSZ)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /STEXP/  BETEXP,SIGMA,ADETLG
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      REAL            GRAD(*),HESS(NP,*),COVAR(NP,*)
      INTEGER         INDX(*)
      NFT2=NFFT/2+1
      CNORM=CBCHI(FITP)
      CALL VRFILL(HESS,0.0,NP*NP)
      CALL VMLTRC(WORK,FR2PIK(1,2),NFT2,FWRK)
      CALL VMLTRC(TWOPIK,FWRK,NFT2,WORK)
      CALL VMLTIC(WORK,NFT2,WORK)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,DDDPAR(1,4))
      CALL FOUR2(WORK,NFFT,1,-1,-1)
      CALL DEGRID(WORK,FWRK)
      CALL VRDOTR(RESID,FWRK,NDAT,HESS(4,4))
      CALL VCOPY(FR2PIK,FWRK,NFFT+2)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,DDDPAR(1,3))
      CALL VCOPY(FR2PIK(1,2),FWRK,NFFT+2)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,WORK)
      CALL VRDOTR(RESID,WORK,NDAT,HESS(3,4))
      HESS(4,3)=HESS(3,4)
      do I=1,NFEW
        J=3+I+I
        AJ=FITP(J)*ASCL
        CALL VMLTRC(EXPF(1,I),FR2PIK,NFT2,WORK)
        CALL VCOPY(WORK,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,J))
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL VMLTB1(FWRK,EXPF(1,I+3),BETEXP,NFT2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,J+1))
        CALL HESS0(HESS,NP,RESID,NDAT,DDDPAR(1,J+1),AJ,J)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL VCOPY(WORK,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,HESS(4,J))
        HESS(J,4)=HESS(4,J)
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL VMLTB1(FWRK,EXPF(1,I+3),BETEXP,NFT2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,SUM)
        HESS(4,J+1)=-AJ*SUM
        HESS(J+1,4)=HESS(4,J+1)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL VMLTB2(WORK,EXPF(1,I+3),BETEXP,NFT2)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SUM)
        HESS(J+1,J+1)=-AJ*SUM
      end do
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC(1,2))
      CALL HESS1(HESS,NP,SCLVEC(1,2),STEPSZ,NFEW)
      CALL MTXINV(HESS,COVAR,NP,INDX,DETLOG)
      END
C     ---------------------------
      SUBROUTINE VMLTB1(F,SK,B,N)
      COMPLEX F(*)
      REAL    SK(*)
      B1=B-1.0
      do I=2,N
        BFACT=B*SK(I)**B1
        F(I)=BFACT*F(I)
      end do
      END
C     ---------------------------
      SUBROUTINE VMLTB2(F,SK,B,N)
      COMPLEX F(*)
      REAL    SK(*)
      DATA    PI /3.141592654/
      PI2=2.0*PI
      B1=B-1.0
      B2=B-2.0
      do I=2,N
        BB1=B*SK(I)**B1
        BFACT=BB1*(BB1-B1/(PI2*SK(I)))
        F(I)=BFACT*F(I)
      end do
      END
