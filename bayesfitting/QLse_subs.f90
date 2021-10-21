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
C     ------------------------------------
      SUBROUTINE OUTPRM(P,C,NP,NFEW,CNORM)
      INCLUDE 'mod_files.f90'
      REAL P(*),C(NP,*)
      IF (NFEW.LT.1 .OR. NFEW.GT.2) RETURN
      OPEN(UNIT=1,FILE=fileout1,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(NFEW,100) P(3),P(1),P(2),P(4)
      do I=5,NP-1,2
        WRITE(NFEW,100) P(I),P(I+1)
      end do
      CSCALE=2.0*CNORM
      do J=1,NP
        do I=1,NP
          C(I,J)=CSCALE*C(I,J)
        end do
      end do
      WRITE(NFEW,100) C(3,3)
      WRITE(NFEW,100) C(3,5),C(5,5)
      WRITE(NFEW,100) C(3,6),C(5,6),C(6,6)
 100  FORMAT(1PE13.4,4E13.4)
      WRITE(NFEW,100) P(NP)
      WRITE(NFEW,100) C(3,NP),C(5,NP),C(6,NP),C(7,NP)
      WRITE(NFEW,*)' -------------------------------------------------'
      close(unit=1)
      END
C
C**<some initialisation routines>***************************************
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
      REAL            V(*)
      PI=3.141592654
      VSMALL=0.0001
      small=0.001
      CCHI=1.0E+10
      CHI=0.0
      B1=BSCL*V(1)
      B2=BSCL*V(2)
      A0=ASCL*V(3)
      DELTAX=V(4)
      NFT2=NFFT/2+1
      CALL CXSHFT(FRES,DELTAX,TWOPIK,FR2PIK,FR2PIK(1,2),NFT2)
      CALL VRFILL(WORK,A0,NFT2)
      IF (NFEW.GT.0) THEN
        BETEXP=V(5+2*NFEW)
        IF (BETEXP.LT.0.01) RETURN
      ENDIF
      PI2=2.0*PI
      do J=1,NFEW
        J3=J+3
        AJ=ASCL*V(3+J+J)
        IF (V(4+J+J).LT.VSMALL) RETURN
        SPNORM=WSCL*V(4+J+J)/(GSCL*PI2)
        EXPF(1,J)=1.0
        EXPF(1,J3)=0.0
        WORK(1,1)=WORK(1,1)+AJ        
        do I=2,NFT2
          SIGIJ=SPNORM*TWOPIK(I)
          EXPIJ=EXP(-PI2*SIGIJ**BETEXP)
          EXPF(I,J)=EXPIJ
          EXPF(I,J3)=SIGIJ
          WORK(I,1)=WORK(I,1)+AJ*EXPIJ
        end do
      end do
      CALL VMLTRC(WORK,FR2PIK,NFT2,FWRK)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,FIT)
      X1=XDAT(1)
      if(o_bgd.eq.2) BNRM=(B2-B1)/(XDAT(NDAT)-X1)
C	avoid conflict BNORM with ModPars
      do I=1,NDAT
        FIT(I)=FIT(I)+B1
         if(o_bgd.eq.2) FIT(I)=FIT(I)+BNRM*(XDAT(I)-X1)
        DIF=FIT(I)-DAT(I)
        RESID(I)=DIF*SIG(I)
        CHI=CHI+DIF*RESID(I)
      end do
      CCHI=CHI/2.0
      sm2=small*small
      v12=V(1)-V(2)
      if(v12.lt.vsmall)v12=1.0
      if(o_bgd.le.1) CCHI=CCHI+v12*v12/sm2
      END
C
C**<search for one more & refine amplitudes>****************************
C
      SUBROUTINE SEARCH(GRAD,HESS,DPAR,NFEW,INDX,COVAR,FITP)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
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
      COMMON /STEXP/  BETEXP,XBETA(101),BETCHI(101),NBTOT
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
        CALL VRDOTR(RESID,WORK(1,2),NDAT,SM)
        HESS(4,J+1)=-AJ*SM
        HESS(J+1,4)=HESS(4,J+1)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL VMLTB2(WORK,EXPF(1,I+3),BETEXP,NFT2)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SM)
        HESS(J+1,J+1)=-AJ*SM
      end do
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC(1,2))
      CALL HESS1(HESS,NP,SCLVEC(1,2),STEPSZ,0)
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
C
C***<see the fit>*******************************************************
C
      SUBROUTINE SEEFIT(SIGPAR,CNORM,PRMSV,SIGSV,NP)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /WRKCOM/ WORK(m_d2,2)
      REAL            SIGPAR(*),PRMSV(*),SIGSV(*)
      CHARACTER       TITLE*80
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      ERRSCL=SQRT(CNORM)
      WRITE(TITLE,100) NFEW
 100  FORMAT(' Best-fit assuming no. of quasi-elastic lines = ',I2)
      WRITE(53,100) NFEW
      WRITE(53,120) CNORM
 120  FORMAT(' >>> Normalised Chi-squared = ',F11.4)
      WRITE(53,130) FITP(1)*BSCL,SIGPAR(1)*BSCL*ERRSCL
 130  FORMAT(' Background(Xmin) = ',1PE12.3,'  +- ',E10.2)
      WRITE(53,140) FITP(2)*BSCL,SIGPAR(2)*BSCL*ERRSCL
 140  FORMAT(' Background(Xmax) = ',1PE12.3,'  +- ',E10.2)
      WRITE(53,150) FITP(4)*GSCL*1000.0,
     1 SIGPAR(4)*GSCL*ERRSCL*1000.0
 150  FORMAT(' Zero offset      = ',F12.2,'  +- ',F10.2,'   ueV')
      WRITE(53,*)' Elastic line'
      WRITE(53,160) FITP(3)*ASCL,SIGPAR(3)*ASCL*ERRSCL
 160  FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
      PRMSV(1)=FITP(3)*ASCL
      SIGSV(1)=SIGPAR(3)*ASCL*ERRSCL
      do I=1,NFEW
        J=4+I+I
        WRITE(53,161) I
 161    FORMAT(' Stretched-exponential : line ',I2)
        WRITE(53,170) 2000.0*FITP(J)*WSCL,
     1   2000.0*SIGPAR(J)*WSCL*ERRSCL
 170    FORMAT(5X,' 2*Sigma    =   ',F13.2,'  +- ',F11.2,'   ueV')
        WRITE(53,180) FITP(J-1)*ASCL,SIGPAR(J-1)*ASCL*ERRSCL
 180    FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
        PRMSV(I+I)=FITP(J-1)*ASCL
        SIGSV(I+I)=SIGPAR(J-1)*ASCL*ERRSCL
        PRMSV(I+I+1)=2.0*FITP(J)*WSCL
        SIGSV(I+I+1)=2.0*SIGPAR(J)*WSCL*ERRSCL
      end do
      IF (NFEW.GT.0) THEN
        BT0=FITP(NP)
        BTSIG=SIGPAR(NP)
        WRITE(53,190) BT0,BTSIG
 190    FORMAT(5X,' Exp-Beta   =   ',F13.3,'  +- ',F11.3)
        PRMSV(NFEW+NFEW+2)=BT0
        SIGSV(NFEW+NFEW+2)=BTSIG
      ENDIF
      close(unit=53)
      END
C     ----------------------------------------------------------
      SUBROUTINE INTDSP(X,Y,XL,XH,YL,YH,SIGPAR,ERRSCL,XMIN,XMAX)
      INCLUDE 'res_par.f90'
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      REAL            X(*),Y(*),XL(*),XH(*),YL(*),YH(*),SIGPAR(*)
      XMAX=0.0
      CALL VRFILL(X,0.0,2*(NFEW+1))
      CALL VRFILL(Y,0.0,2*(NFEW+1))
      Y(2)=FITP(3)*ASCL
      XL(1)=0.0
      XH(1)=0.0
      YL(1)=Y(2)-ASCL*SIGPAR(3)*ERRSCL
      YH(1)=Y(2)+ASCL*SIGPAR(3)*ERRSCL
      do I=1,NFEW
        J=I+I+1
        X(J)=2.0*FITP(J+3)*WSCL
        X(J+1)=X(J)
        Y(J+1)=FITP(J+2)*ASCL
        SIGX=SIGPAR(J+3)*WSCL*ERRSCL*2.0
        SIGY=SIGPAR(J+2)*ASCL*ERRSCL
        XL(I+1)=X(J+1)-SIGX
        XH(I+1)=X(J+1)+SIGX
        YL(I+1)=Y(J+1)-SIGY
        YH(I+1)=Y(J+1)+SIGY
        IF (XH(I+1).GT.XMAX) XMAX=XH(I+1)
      end do
      IF (NFEW.EQ.0) XMAX=WSCL/10.0
      XMIN=-XMAX/50.0
      XMAX=XMAX*1.03
      END
C
