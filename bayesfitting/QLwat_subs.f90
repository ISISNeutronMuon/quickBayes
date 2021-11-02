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
      END
C     ------------------------------------
      SUBROUTINE OUTPRM(P,C,NP,NFEW,CNORM)
      INCLUDE 'mod_files.f90'
      REAL P(*),C(NP,*)
      IF (NFEW.NE.1) RETURN
       OPEN(UNIT=1,FILE=fileout1,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(NFEW,100) P(3),P(1),P(2),P(4)
      WRITE(NFEW,100) P(5),P(6),P(7)
      CSCALE=2.0*CNORM
      do J=1,NP
       do I=1,NP
        C(I,J)=CSCALE*C(I,J)
       end do
      end do
      WRITE(NFEW,100) C(3,3)
      WRITE(NFEW,100) C(3,5),C(5,5)
      WRITE(NFEW,100) C(3,6),C(5,6),C(6,6)
      WRITE(NFEW,100) C(3,7),C(5,7),C(6,7),C(7,7)
 100  FORMAT(1PE13.4,4E13.4)
      WRITE(NFEW,*)' -------------------------------------------------'
      close(unit=1)
      return
      END
c
C***<calculate data & chi-squared>**************************************
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
      COMMON /BESSJ/  AJ0,AJ1,AJ2
      REAL   V(*)
      DATA   E3,E9 /3.0,9.0/
      PI=3.141592654
      VSMALL=0.0001
      small=0.001
      CCHI=1.0E+10
      CHI=0.0
      BNRM=1.0
      B1=BSCL*V(1)
      B2=BSCL*V(2)
      A0=ASCL*V(3)
      DELTAX=V(4)
      NFT2=NFFT/2+1
      CALL CXSHFT(FRES,DELTAX,TWOPIK,FR2PIK,FR2PIK(1,2),NFT2)
      CALL VRFILL(WORK,A0,NFT2)
      XNSCL=-LOG(1.0E-7)/(TWOPIK(2)-TWOPIK(1))
      IF (NFEW.EQ.1) THEN
        E3B=1.0/E3
        E9B=1.0/E9
        A1=ASCL*V(5)
        SIG1=WSCL*V(6)/GSCL
        SIG2=WSCL*V(7)/(GSCL*3.0)
        IF (SIG1.LT.0.0 .OR. SIG2.LT.0.0) RETURN
        CALL VRFILL(EXPF(1,1),0.0,NFT2)
        CALL VRFILL(EXPF(1,2),0.0,NFT2)
        CALL VRFILL(EXPF(1,3),0.0,NFT2)
        NXFT21=1+NINT(XNSCL/(ABS(SIG1)+1.0E-10))
        NXFT22=1+NINT(XNSCL/(ABS(SIG2)+1.0E-10))
        NXFT2=MIN(MAX(NXFT21,NXFT22),NFT2)
        do I=1,NXFT2
          EXPIJ1=EXP(-TWOPIK(I)*SIG1)
          EXPIJ2=EXP(-TWOPIK(I)*SIG2)
          EXPIJ=EXPIJ1*(AJ0+EXPIJ2*(AJ1+AJ2*EXPIJ2**2))
          WORK(I,1)=WORK(I,1)+A1*EXPIJ
          EXPF(I,1)=EXPIJ
          EXPF(I,2)=E3B*EXPIJ1*EXPIJ2*(AJ1+E3*AJ2*EXPIJ2**2)
          EXPF(I,3)=E9B*EXPIJ1*EXPIJ2*(AJ1+E9*AJ2*EXPIJ2**2)
        end do
      ENDIF
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
      CCHI=CHI/(2.0*FLOAT(NDAT))
      sm2=small*small
      v12=V(1)-V(2)
      if(v12.lt.vsmall)v12=1.0
      if(o_bgd.le.1) CCHI=CCHI+v12*v12/sm2
      RETURN
      END
C
C**<search for one more & refine amplitudes>****************************
C
      SUBROUTINE SEARCHw(GRAD,HESS,DPAR,NFEW,INDX,COVAR,FITP)
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      REAL     GRAD(*),HESS(*),DPAR(*),COVAR(*),FITP(*)
      INTEGER  INDX(*)
      J=4+3*NFEW
      DXLOG=0.85
      NSRCH=NINT(LOG(5.0*GSCL/WSCL)/LOG(DXLOG))
      CMIN=1.0E20
      SIG1=0.
      SIG2=0.
      FITP(J-1)=1.0
      do I2=1,NSRCH
        FITP(J)=1.0
        FITP(J-2)=0.5
        do I1=1,NSRCH
          CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
          CNORM=CCHI(FITP)
          IF (CNORM.LT.CMIN) THEN
            CMIN=CNORM
            SIG1=FITP(J-1)
            SIG2=FITP(J)
          ENDIF
          FITP(J)=FITP(J)*DXLOG
        end do
        FITP(J-1)=FITP(J-1)*DXLOG
      end do
      FITP(J-1)=SIG1
      FITP(J)=SIG2
      CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
      END
C**<refinement>*********************************************************
      SUBROUTINE REFINAw(GRAD,HESS,DPAR,NP,DETLOG,INDX,COVAR)
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
      CNORM=CCHI(FITP)
      CALL VRFILL(HESS,0.0,NP*NP)
      CALL VCOPY(FR2PIK,FWRK,NFFT+2)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,DDDPAR(1,3))
      if (NFEW.eq.1)then
        CALL VMLTRC(EXPF(1,1),FR2PIK,NFT2,FWRK)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,4))
      endif
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC)
      CALL HESS1(HESS,NP,SCLVEC,0.3,NFEW)
      CALL INVERT(HESS,COVAR,NP,INDX,DETLOG)
      CALL NEWEST(COVAR,GRAD,NP,NFEW,DPAR,FITP)
      CNORM=CCHI(FITP)
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC)
      CALL NEWEST(COVAR,GRAD,NP,NFEW,DPAR,FITP)
      END
C     -------------------------------------------------------
      SUBROUTINE REFINEw(GRAD,HESS,NP,DETLOG,INDX,COVAR,STEPSZ) 
      INCLUDE 'res_par.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      COMMON /BESSJ/  AJ0,AJ1,AJ2
      REAL GRAD(*),HESS(NP,*),COVAR(NP,*)
      INTEGER INDX(*)
C
      NFT2=NFFT/2+1
      CNORM=CCHI(FITP)
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
      IF (NFEW.EQ.1) THEN
        A1=FITP(5)*ASCL
        CALL VMLTRC(EXPF(1,1),FR2PIK,NFT2,WORK)
        CALL VCOPY(WORK,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,5))
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,6))
        CALL HESS0(HESS,NP,RESID,NDAT,DDDPAR(1,6),A1,5)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL VCOPY(WORK,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,HESS(4,5))
        HESS(5,4)=HESS(4,5)
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,SM)
        HESS(4,6)=-A1*SM
        HESS(6,4)=HESS(4,6)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SM)
        HESS(6,6)=-A1*SM
        CALL VMLTRC(EXPF(1,2),FR2PIK,NFT2,WORK)
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,7))
        CALL HESS0(HESS,NP,RESID,NDAT,DDDPAR(1,7),A1,5)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,SM)
        HESS(4,7)=-A1*SM
        HESS(7,4)=HESS(4,7)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SM)
        HESS(6,7)=-A1*SM
        HESS(7,6)=HESS(6,7)
        CALL VMLTRC(EXPF(1,3),FR2PIK,NFT2,WORK)
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VMLTRC(TWOPIK,FWRK,NFT2,WORK)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SM)
        HESS(7,7)=A1*SM
      ENDIF
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC(1,2))
      CALL HESS1(HESS,NP,SCLVEC(1,2),STEPSZ,0)
      CALL INVERT(HESS,COVAR,NP,INDX,DETLOG)
      RETURN
      END
C***<see the fit>*******************************************************
      SUBROUTINE SEEFITw(SIGPAR,CNORM,PRMSV,SIGSV)
      INCLUDE 'res_par.f90'
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /WRKCOM/ WORK(m_d2,2)
      INCLUDE 'mod_files.f90'
      REAL       SIGPAR(*),PRMSV(*),SIGSV(*)
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      ERRSCL=SQRT(CNORM)
      NQUASI=0
      IF (NFEW.EQ.1) NQUASI=3
      WRITE(53,100) NQUASI
 100  FORMAT(' Best-fit assuming no. of quasi-elastic lines = ',I2)
      WRITE(53,120) CNORM
 120  FORMAT(' Normalised Chi-squared = ',F11.4)
      WRITE(53,130) FITP(1)*BSCL,SIGPAR(1)*BSCL*ERRSCL
 130  FORMAT(' Background(Xmin) = ',1PE12.3,'  +- ',E10.2)
      WRITE(53,99) SIGPAR(1),BSCL,ERRSCL
  99  FORMAT(' Sigpar = ',3(1PE12.3,1x))
      WRITE(53,140) FITP(2)*BSCL,SIGPAR(2)*BSCL*ERRSCL
 140  FORMAT(' Background(Xmax) = ',1PE12.3,'  +- ',E10.2)
      WRITE(53,150) FITP(4)*GSCL*1000.0,
     1  SIGPAR(4)*GSCL*ERRSCL*1000.0
 150  FORMAT(' Zero offset       = ',F12.2,'  +- ',F10.2,'  ueV')
      WRITE(53,*)' Elastic line' 
      WRITE(53,160) FITP(3)*ASCL,SIGPAR(3)*ASCL*ERRSCL
 160  FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
      PRMSV(1)=FITP(3)*ASCL
      SIGSV(1)=SIGPAR(3)*ASCL*ERRSCL
      IF (NFEW.EQ.1) THEN
        WRITE(53,161)
 161    FORMAT(' Quasi-elastic parameters')
        WRITE(53,170) 2.0*FITP(6)*WSCL,2.0*SIGPAR(6)*WSCL*ERRSCL
 170    FORMAT(5X,' FWHM0      =   ',F13.2,'  +- ',F11.2)
        WRITE(53,171) 2.0*FITP(7)*WSCL,2.0*SIGPAR(7)*WSCL*ERRSCL
 171    FORMAT(5X,' FWHM1      =   ',F13.2,'  +- ',F11.2)
        WRITE(53,180) FITP(5)*ASCL,SIGPAR(5)*ASCL*ERRSCL
 180    FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
        PRMSV(2)=FITP(5)*ASCL
        SIGSV(2)=SIGPAR(5)*ASCL*ERRSCL
        PRMSV(3)=2.0*FITP(6)*WSCL
        SIGSV(3)=2.0*SIGPAR(6)*WSCL*ERRSCL
        PRMSV(4)=2.0*FITP(7)*WSCL
        SIGSV(4)=2.0*SIGPAR(7)*WSCL*ERRSCL
      ENDIF
      close(unit=53)
      END
C     ----------------------------------
      SUBROUTINE ERRBARw(COVAR,NP,SIGPAR,COV12)
      REAL COVAR(NP,*),SIGPAR(*)
      SMALL=1.0E-20
      do I=1,NP
        SIGPAR(I)=SQRT(2.0*ABS(COVAR(I,I))+SMALL)
      end do
      IF (NP.EQ.7) COV12=2.0*COVAR(6,7)
      END
C     ----------------------------------------------------------
      SUBROUTINE INTDSPw(X,Y,XL,XH,YL,YH,SIGPAR,ERRSCL,XMIN,XMAX,COV12)
      INCLUDE 'res_par.f90'
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /BESSJ/  AJ0,AJ1,AJ2
      REAL   X(*),Y(*),XL(*),XH(*),YL(*),YH(*),SIGPAR(*)
      XMAX=0.0
      CALL VRFILL(X,0.0,2*(NFEW+1))
      CALL VRFILL(Y,0.0,2*(NFEW+1))
      Y(2)=FITP(3)*ASCL
      XL(1)=0.0
      XH(1)=0.0
      YL(1)=Y(2)-ASCL*SIGPAR(3)*ERRSCL
      YH(1)=Y(2)+ASCL*SIGPAR(3)*ERRSCL
      IF (NFEW.EQ.1) THEN
       A1=FITP(5)*ASCL
       FWHM1=FITP(6)*WSCL*2.0
       FWHM2=FITP(7)*WSCL*2.0
       SIGA1=SIGPAR(5)*ASCL*ERRSCL
       SIGFW1=SIGPAR(6)*WSCL*ERRSCL*2.0
       SIGFW2=SIGPAR(7)*WSCL*ERRSCL*2.0
       COV12B=COV12*(WSCL*ERRSCL*2.0)**2
       CALL VRFILL(X(3),FWHM1,2)
       Y(4)=A1*AJ0
       SIGX=SIGFW1
       SIGY=SIGA1*AJ0
       XL(2)=X(4)-SIGX
       XH(2)=X(4)+SIGX
       YL(2)=Y(4)-SIGY
       YH(2)=Y(4)+SIGY
       CALL VRFILL(X(5),FWHM1+FWHM2/3.0,2)
       Y(6)=A1*AJ1
       SIGX=SQRT(SIGFW1**2+(SIGFW2**2)/9.0+2.0*COV12B/3.0)
       SIGY=SIGA1*AJ1
       XL(3)=X(6)-SIGX
       XH(3)=X(6)+SIGX
       YL(3)=Y(6)-SIGY
       YH(3)=Y(6)+SIGY
       CALL VRFILL(X(7),FWHM1+FWHM2,2)
       Y(8)=A1*AJ2
       SIGX=SQRT(SIGFW1**2+SIGFW2**2+2.0*COV12B)
       SIGY=SIGA1*AJ2
       XL(4)=X(8)-SIGX
       XH(4)=X(8)+SIGX
       YL(4)=Y(8)-SIGY
       YH(4)=Y(8)+SIGY
       XMAX=MAX(XMAX,XH(2),XH(3),XH(4))
       XMIN=-XMAX/50.0
       XMAX=XMAX*1.03
      ELSE
       XMAX=WSCL/20.0
       XMIN=-XMAX/5.0
      ENDIF
      END
C     ------------------------------------------------------------------
      SUBROUTINE PROBNw(CNORM,NDAT,DETLOG,NFEW,PROBLG)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      EXTERNAL FCTNLG
      CHISQ=CNORM*FLOAT(NDAT)
      DETLOG=DETLOG-FLOAT(NFEW*2)*LOG10(ASCL*WSCL)
      CHI2LG=-LOG10(2.7182818)*CHISQ/2.0
      PROBLG=CHI2LG-(0.5*DETLOG)-(FLOAT(NFEW)*LOG10(ASCL*WSCL))
      PROBLG=PROBLG+FCTNLG(NFEW)
      PROBLG=PROBLG+(FLOAT(NFEW)*LOG10(4.0*3.141592654))
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(53,110) PROBLG
 110  FORMAT(' Log10[Prob(Stretched exp|{Data})] = ',F11.1)
      IF (NFEW.LT.1) WRITE(53,*)' -------------------------'
      close(unit=53)
      END
C     ---------------------------------
      SUBROUTINE PRBOUTw(P,NP,NQ,POUT)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      REAL  P(4,m_sp),POUT(4,m_sp)
      J=NQ
      SM=0.0
      do I=1,NP
        SM=SM+10.0**(P(I,J)-P(NP,J))
      end do
      PLNORM=LOG10(SM)+P(NP,J)
      do I=1,NP
        P(I,J)=P(I,J)-PLNORM
      end do
      do I=1,NP
        POUT(I,J)=P(I,J)
       end do
      END
C     ---------------------------------
      SUBROUTINE BESINT(AJ0,AJ1,AJ2,QA)
      INCLUDE 'mod_files.f90'
      BJ0=BESSJ0(QA)
      BJ1=BESSJ1(QA)
      BJ2=BJ1*2.0/QA-BJ0
      AJ0=BJ0**2
      AJ1=3.0*BJ1**2
      AJ2=5.0*BJ2**2
      RETURN
      END
C     ---------------------------------
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
      RETURN
      END
C***<Numerical Recipes>*************************************************
      FUNCTION bessj0(x)
      REAL bessj0,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *                 s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *     s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     * -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     * .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     * 651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     * s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     * 59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *         (s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *              p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      END
C     ------------------
      FUNCTION bessj1(x)
      REAL bessj1,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *                 s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *     s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     * 242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     * s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     * 99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     * .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     * -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *         y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *    p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      endif
      END
C
