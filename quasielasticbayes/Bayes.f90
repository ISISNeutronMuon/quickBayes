C***<read in the data>**************************************************
      SUBROUTINE DATIN(IREAD,DTNORM,IDUF)
      INCLUDE 'mod_files.f90'
      INCLUDE 'mod_data.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      REAL DTNORM(m_sp)
      IDUF=0
      SMALL=1.0E-10
      DSUM=0.0
      COS2TH=2.0*COS(theta(IREAD)*3.14159265/180.0)
      QQ=efix+efix-COS2TH*SQRT(efix*efix)
      QAVRG(IREAD)=0.69469*SQRT(QQ)
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(53,*)' ----------------------------------------------------'
      write(53,1001)IREAD,theta(IREAD),QAVRG(IREAD)
1001  format(' Group ',i2,' : theta = ',f10.5,' Q = ',f10.5)
      WRITE(53,*)' ----------------------------------------------------'
      DTNRM=DTNORM(IREAD)
      NDAT=ntc-1
      do I=1,NDAT
       XDAT(I)=xin(I)
       DAT(I)=yin(I)*DTNRM
       SIG(I)=ein(I)*DTNRM
       IF (SIG(I).GT.SMALL) THEN
        SIG(I)=SIG(I)*SIG(I)
        DSUM=DSUM+DAT(I)
       ELSE
        SIG(I)=0.0
       ENDIF
      end do
      IF (DSUM.LT.SMALL) IDUF=1
      CALL DATIN1
      close(unit=53)
      END
C     ---------------------------------------------------------
      SUBROUTINE DATIN1
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON/ModPars/NBIN,IMIN,IMAX,RSCL,BNORM
      SMALL=1.0E-20
      if (ABS(RSCL-1.0).GT.0.01)then
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       WRITE(53,*) ' DATIN1; Data error-bars multiplied by: ',RSCL
       close(unit=53)
      endif
      RSCL=RSCL*RSCL
      N=0
      do I=IMIN,IMAX,NBIN
       N=N+1
       XXD=0.0
       DD=0.0
       EE=0.0
       K=0
       do J=0,NBIN-1
        XXD=XXD+XDAT(I+J)
        IF (SIG(I+J).GT.SMALL) THEN
         K=K+1
         DD=DD+DAT(I+J)
         EE=EE+SIG(I+J)
        ENDIF
       end do
       XDAT(N)=BNORM*XXD
       IF (K.GT.0) THEN
        DAT(N)=BNORM*DD
        SIG(N)=2.0*FLOAT(K*K)/(EE*RSCL)
       ELSE
        DAT(N)=0.0
        SIG(N)=0.0
       ENDIF
      end do
      NDAT=N
      END
C     ------------------
      FUNCTION FCTNLG(N)
      A=0.0
      do I=1,N
       A=A+LOG10(FLOAT(I))
      end do
      FCTNLG=A
      END
C     ----------------------------------------------
      SUBROUTINE CXSHFT(RK,DX,TWOPIK,RKEXP,RKEXP2,N)
      COMPLEX RK(*),RKEXP(*),RKEXP2(*),XC
      REAL    TWOPIK(*)
      do I=1,N
        XX=TWOPIK(I)*DX
        XC=CMPLX(COS(XX),SIN(XX))
        RKEXP(I)=RK(I)*XC
      end do
      CALL VMLTRC(TWOPIK,RKEXP,N,RKEXP2)
      CALL VMLTIC(RKEXP2,N,RKEXP2)
      END
C     ---------------------------
      SUBROUTINE VMLTRC(R,C,N,CC)
      REAL    R(*)
      COMPLEX C(*),CC(*)
      do I=1,N
        A=R(I)*REAL(C(I))
        B=R(I)*AIMAG(C(I))
        CC(I)=CMPLX(A,B)
      end do
      END
C     -------------------------
      SUBROUTINE VMLTIC(C,N,CI)
      COMPLEX C(*),CI(*)
      do J=1,N
        A=REAL(C(J))
        B=AIMAG(C(J))
        CI(J)=CMPLX(-B,A)
      end do
      END
C     ---------------------------
      SUBROUTINE MLTMXV(P,OP,N,D)
      REAL P(*),OP(N,N),D(*)
      do K=1,N
        SM=0.0
        do J=1,N
          SM=SM+OP(J,K)*P(J)
        end do
        D(K)=SM
      end do
      END
C     ---------------------------------------------
      SUBROUTINE DEGRID(YGRD,YDAT)
      INCLUDE 'res_par.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      REAL    YGRD(*),YDAT(*)
      do I=1,NDAT
        J=IPDAT(I)
        YDAT(I)=YGRD(J)+XPDAT(I)*(YGRD(J+1)-YGRD(J))
      end do
      END
C     --------------------------------------------
      SUBROUTINE GRADPR(GRAD,RESID,NDAT,NP,SCLVEC)
      INCLUDE 'res_par.f90'
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      REAL   GRAD(*),RESID(*),SCLVEC(*)
      do I=1,NP
        CALL VRDOTR(RESID,DDDPAR(1,I),NDAT,SM)
        GRAD(I)=SCLVEC(I)*SM
      end do
      END
C     --------------------------
      SUBROUTINE VRDOTR(A,B,N,C)
      REAL A(*),B(*)
      C=0.0
      do I=1,N
       C=C+A(I)*B(I)
      end do
      END
C     ------------------------------------------------
      SUBROUTINE HESS0(HESS,NP,RESID,NDAT,DDDPAR,AJ,J)
      REAL HESS(NP,*),RESID(*),DDDPAR(*)
      SM=0.0
      do I=1,NDAT
        SM=SM-RESID(I)*DDDPAR(I)
        DDDPAR(I)=-AJ*DDDPAR(I)
      end do
      HESS(J,J+1)=SM
      HESS(J+1,J)=SM
      END
C     ----------------------------------
      SUBROUTINE ERRBAR(COVAR,NP,SIGPAR)
      REAL COVAR(NP,*),SIGPAR(*)
      SMALL=1.0E-20
      do I=1,NP
        SIGPAR(I)=SQRT(2.0*ABS(COVAR(I,I))+SMALL)
      end do
      END
C     ----------------------------------------
      SUBROUTINE DPINIT
      INCLUDE 'res_par.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      I=1
      GSCL=(XJ(NFFT)-XJ(1))/FLOAT(NFFT-1)
      do K=1,NDAT
        do J=I,NFFT
          IF (XJ(J).GE.XDAT(K)) GOTO 1
        end do
   1    IPDAT(K)=J-1
        XPDAT(K)=(XDAT(K)-XJ(J-1))/GSCL
        I=J
      end do
      END
C     -----------------
      SUBROUTINE GDINIT
      INCLUDE 'res_par.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      X1=XDAT(1)
      XN=XDAT(NDAT)
      GNORM=1.0/(XN-X1)
      do I=1,NDAT
        XGNORM=(XDAT(I)-X1)*GNORM
        DDDPAR(I,1)=1.0-XGNORM
        DDDPAR(I,2)=XGNORM
      end do
      END
C     ---------------------------------------------
      SUBROUTINE INVERT(HESS,COVAR,NP,INDX,DETLOG)
      INCLUDE 'options.f90'
      REAL    HESS(NP,*),COVAR(NP,*)
      INTEGER INDX(*)
      SMALL=1.E-20
      DETLOG=0.0
      if (prog.eq.'s') then
       cov=2.0
      else
       cov=1.0
      endif
      CALL VRFILL(COVAR,0.0,NP*NP)
      do I=1,NP
       COVAR(I,I)=cov
      end do
      CALL LUDCMP(HESS,NP,NP,INDX,D)
      do I=1,NP
        DETLOG=DETLOG+LOG10(ABS(HESS(I,I))+SMALL)
        IF (HESS(I,I).LT.0.0) D=-D
      end do
      do I=1,NP
      CALL LUBKSB(HESS,NP,NP,INDX,COVAR(1,I))
      end do
      END
C     ---------------------------------
      SUBROUTINE PRBOUT(P,NP,NQ,POUT)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      REAL  P(4,m_sp),POUT(4,m_sp)
      J=NQ
      SM=0.0
      do I=1,NP
        SM=SM+10.0**(P(I,J)-P(NP,J))
      end do
      PLMAX=MAX(P(3,J),P(4,J))
      PLNORM=LOG10(SM)+P(NP,J)
      do I=1,NP
        P(I,J)=P(I,J)-PLMAX
      end do
      do I=1,NP
        POUT(I,J)=P(I,J)
       end do
      END
C     -----------------------------------------------
      SUBROUTINE NEWEST(COVAR,GRAD,NP,NFEW,DPAR,FITP)
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      REAL COVAR(NP,*),GRAD(*),DPAR(*),FITP(*)
      if (prog.eq.'w') then
       mp=3
      else
       mp=2
      endif
      CALL MLTMXV(GRAD,COVAR,NP,DPAR)
      IF (NP.EQ.(4+mp*NFEW)) THEN
        do I=1,NP
          FITP(I)=FITP(I)-DPAR(I)
        end do
      ELSEIF(NP.EQ.(3+NFEW)) THEN
        do I=1,3
          FITP(I)=FITP(I)-DPAR(I)
        end do
        do I=1,NFEW
          J=I+3
          FITP(J+I)=FITP(J+I)-DPAR(J)
        end do
      ELSE
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       WRITE(53,*)' NEWEST : Something wrong here folks!'
       STOP 
       close(unit=53)
      ENDIF
      END
C     ------------------------------------------------
      SUBROUTINE HESS1(HESS,NP,SCLVEC,STEPSZ,NFEW)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      REAL   HESS(NP,*),SCLVEC(*)
      do J=1,NP
        do I=J,NP
          SM=0.0
          do K=1,NDAT
           SM=SM+SIG(K)*DDDPAR(K,I)*DDDPAR(K,J)
          end do
          HESS(I,J)=(HESS(I,J)+SM)*SCLVEC(I)*SCLVEC(J)
          HESS(J,I)=HESS(I,J)
        end do
      end do
      BEEFUP=2.0/(STEPSZ*STEPSZ)
      do I=1,NP
        HESS(I,I)=HESS(I,I)+BEEFUP
      end do
      if (prog.eq.'l'.OR.prog.eq.'s') then
       IF (NFEW.GT.0) then                        !option for elastic peak
        if (o_el.eq.0) HESS(3,3)=2.0E8
       ENDIF
      endif
      END
C     ------------------------------------------------------------------
      SUBROUTINE PROBN(CNORM,NDAT,DETLOG,NFEW,NMAX,PRBSV)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      REAL     PRBSV(*)
      EXTERNAL FCTNLG
      CHISQ=CNORM*FLOAT(NDAT)
      DETLOG=DETLOG-FLOAT(NFEW*2)*LOG10(ASCL*WSCL)
      CHI2LG=-LOG10(2.7182818)*CHISQ/2.0
      PROBLG=CHI2LG-(0.5*DETLOG)-(FLOAT(NFEW)*LOG10(ASCL*WSCL))
      PROBLG=PROBLG+FCTNLG(NFEW)
      PROBLG=PROBLG+(FLOAT(NFEW)*LOG10(4.0*3.141592654))
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      if (prog.eq.'l') then
       WRITE(53,101) NFEW,PROBLG
 101   FORMAT(' Log10[Prob(',I2,' Quasi-elastic lines|{Data})] = ',
     1  F11.1)
       IF (NFEW.LT.NMAX) WRITE(53,*)' -------------------------'
      endif
      if (prog.eq.'s') then
       WRITE(53,102)  PROBLG
 102   FORMAT(' Log10[Prob(Stretched exp|{Data})] = ',F11.1)
       IF (NFEW.LT.1) WRITE(53,*)' -------------------------'
      endif
      if (prog.eq.'w') then
       WRITE(53,103)  PROBLG
 103   FORMAT(' Log10[Prob(Water|{Data})] = ',F11.1)
       IF (NFEW.LT.1) WRITE(53,*)' -------------------------'
      endif
      close(unit=53)
      PRBSV(NFEW+1)=PROBLG
      END
C     -------------------------------------------------------
      SUBROUTINE REFINA(GRAD,HESS,DPAR,NP,DETLOG,INDX,COVAR)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      REAL            GRAD(*),HESS(NP,*),DPAR(*),COVAR(NP,*)
      INTEGER         INDX(*)
      NFT2=NFFT/2+1
      if (prog.eq.'s'.OR.prog.eq.'q') then
       CNORM=CBCHI(FITP)
      else
       CNORM=CCHI(FITP)
      endif
      CALL VRFILL(HESS,0.0,NP*NP)
      CALL VCOPY(FR2PIK,FWRK,NFFT+2)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,DDDPAR(1,3))
      do I=1,NFEW
        CALL VMLTRC(EXPF(1,I),FR2PIK,NFT2,FWRK)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,3+I))
      end do
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC)
      CALL HESS1(HESS,NP,SCLVEC,0.3,NFEW)
      if (prog.eq.'s'.OR.prog.eq.'q') then
       CALL MTXINV(HESS,COVAR,NP,INDX,DETLOG)
      else
       CALL INVERT(HESS,COVAR,NP,INDX,DETLOG)
      endif
      CALL NEWEST(COVAR,GRAD,NP,NFEW,DPAR,FITP)
      if (prog.eq.'s'.OR.prog.eq.'q') then
       CNORM=CBCHI(FITP)
      else
       CNORM=CCHI(FITP)
      endif
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC)
      CALL NEWEST(COVAR,GRAD,NP,NFEW,DPAR,FITP)
      END
C     ---------------------------------------------
      SUBROUTINE MTXINV(HESS,COVAR,NP,INDX,DETLOG)
      REAL    HESS(NP,*),COVAR(NP,*)
      INTEGER INDX(*)
      SMALL=1.E-20
      DETLOG=0.0
      CALL VRFILL(COVAR,0.0,NP*NP)
      do I=1,NP
        COVAR(I,I)=1.0
      end do
      CALL LUDCMP(HESS,NP,NP,INDX,D)
      do I=1,NP
        DETLOG=DETLOG+LOG10(ABS(HESS(I,I))+SMALL)
        IF (HESS(I,I).LT.0.0) D=-D
      end do
      do I=1,NP
        CALL LUBKSB(HESS,NP,NP,INDX,COVAR(1,I))
      end do
      END
C     ---------------------------------------------
      FUNCTION CBCHI(V)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /STEXP/  BETEXP,XBETA(101),BETCHI(101),NBTOT
      REAL   V(*)
      PI=3.141592654
      VSMALL=0.0001
      CBCHI=1.0E+10
      CHI=0.0
      B1=BSCL*V(1)
      B2=BSCL*V(2)
      A0=ASCL*V(3)
      DELTAX=V(4)
      NFT2=NFFT/2+1
      CALL CXSHFT(FRES,DELTAX,TWOPIK,FR2PIK,FR2PIK(1,2),NFT2)
      CALL VRFILL(WORK,A0,NFT2)
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
C     avoid conflict BNORM with ModPars
      do I=1,NDAT
        FIT(I)=FIT(I)+B1
         if(o_bgd.eq.2) FIT(I)=FIT(I)+BNRM*(XDAT(I)-X1)
        DIF=FIT(I)-DAT(I)
        RESID(I)=DIF*SIG(I)
        CHI=CHI+DIF*RESID(I)
      end do
      CBCHI=CHI/(2.0*FLOAT(NDAT))
      END
C     ---------------------------------------------
      SUBROUTINE PRINIT(FITP,NQMAX,NFEW,IXSCAL)
      INCLUDE 'mod_files.f90'
      INCLUDE 'mod_data.f90'
      INCLUDE 'options.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      REAL FITP(m_p)
      if (IXSCAL.gt.1) goto 3
      SMALL=1.0E-10
      SM=0.0
      NSUM=0
      do I=1,20
       IF (SIG(I).GE.SMALL) THEN
        NSUM=NSUM+1
        SM=SM+ABS(DAT(I))
       ENDIF
       IF (NSUM.GE.10) GOTO 1
      end do
   1  BSCL1=SM/FLOAT(NSUM)
      SM=0.0
      NSUM=0
      do I=1,20
       IF (SIG(NDAT-I+1).GE.SMALL) THEN
        NSUM=NSUM+1
        SM=SM+ABS(DAT(NDAT-I+1))
       ENDIF
       IF (NSUM.GE.10) GOTO 2
      end do
   2  BSCL=SM/FLOAT(NSUM)
      IF (BSCL1.LT.BSCL) BSCL=BSCL1
      BSCL=BSCL/2.0
      if(o_bgd.eq.0)BSCL=0.0            !zero background
      AXMAX=ABS(XDAT(NDAT))
      IF (ABS(XDAT(1)).GT.AXMAX) AXMAX=ABS(XDAT(1))
      WSCL=AXMAX/3.0
      MK=0
      SM=0.0
      SUMSIG=0.0
      do I=1,NDAT-1
       IF (SIG(I).GE.SMALL) THEN
        MK=MK+1
        SM=SM+(DAT(I)-BSCL)*(XDAT(I+1)-XDAT(I))
        SUMSIG=SUMSIG+SQRT((2.0/SIG(I)))
       ENDIF
      end do
      ASCL=SM*FLOAT(NDAT)/FLOAT(MK)
      SUMSIG=SUMSIG/FLOAT(MK)
      SUMSIG=(XDAT(NDAT)-XDAT(1))*SUMSIG/SQRT(FLOAT(MK))
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      IF (ASCL.LT.SUMSIG) THEN
       write(53,1001)
1001   format(' qlm> *** Estimate of Amax is being set to lower bound!')
       write(53,1002)ASCL,SUMSIG
1002   format(' ( ',e14.7,' --> ',e14.7,' )')
       ASCL=SUMSIG
      ENDIF
      WRITE(53,*)' ----------------------------------------------------'
      close(unit=53)
      ASCL=ASCL/GSCL
      FITP(1)=1.0
      FITP(2)=1.0
      FITP(3)=0.5
      do I=1,2
       SCLVEC(I,1)=BSCL
       SCLVEC(I,2)=BSCL
      end do
      SCLVEC(3,1)=ASCL
      SCLVEC(3,2)=ASCL
    3 SCLVEC(4,2)=1.0
      if(prog.eq.'w')then
       SCLVEC(4,1)=ASCL
       SCLVEC(5,2)=ASCL
       SCLVEC(6,2)=WSCL/GSCL
       SCLVEC(7,2)=WSCL/GSCL
      else
       do I=1,NQMAX
        SCLVEC(3+I,1)=ASCL
        SCLVEC(3+I+I,2)=ASCL
        SCLVEC(4+I+I,2)=WSCL/GSCL
       end do
      endif
      NFEW=0
      END
C     ---------------------------------------------
      SUBROUTINE FileInit(Nu,ISP)
      INCLUDE 'mod_files.f90'
      INCLUDE 'mod_data.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      OPEN(UNIT=1,FILE=fileout1,STATUS='old',FORM='formatted',
     1 access='append')
      if (Nu.ge.2) then
       OPEN(UNIT=2,FILE=fileout2,STATUS='old',FORM='formatted',
     1 access='append')
       if (Nu.eq.3) then
        OPEN(UNIT=3,FILE=fileout3,STATUS='old',FORM='formatted',
     1  access='append')
       endif
      endif
      do n=1,Nu
       write(n,*) QAVRG(ISP),ASCL,WSCL,BSCL,GSCL
       close(unit=n)
      end do
      END

