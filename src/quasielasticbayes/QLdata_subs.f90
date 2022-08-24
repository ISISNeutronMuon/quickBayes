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
      IF (NFEW.LT.1 .OR. NFEW.GT.3) RETURN
      OPEN(UNIT=1,FILE=fileout1,STATUS='old',FORM='formatted',
     1 access='append')
      OPEN(UNIT=2,FILE=fileout2,STATUS='old',FORM='formatted',
     1 access='append')
      OPEN(UNIT=3,FILE=fileout3,STATUS='old',FORM='formatted',
     1 access='append')
      WRITE(NFEW,100) P(3),P(1),P(2),P(4)
      do I=5,NP,2
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
      IF (NFEW.GT.1) THEN
        WRITE(NFEW,100) C(3,7),C(5,7),C(6,7),C(7,7)
        WRITE(NFEW,100) C(3,8),C(5,8),C(6,8),C(7,8),C(8,8)
      ENDIF
      IF (NFEW.GT.2) THEN
        WRITE(NFEW,100) C(3,9),C(5,9),C(6,9),C(7,9),C(8,9),C(9,9)
        WRITE(NFEW,100) C(3,10),C(5,10),C(6,10),C(7,10),C(8,10),
     1 C(9,10),C(10,10)
      ENDIF
 100  FORMAT(1PE13.4,6E13.4)
      WRITE(NFEW,*)' -------------------------------------------------'
      close(unit=1)
      close(unit=2)
      close(unit=3)
      END
C
C***<set up blur function>**********************************************

      SUBROUTINE BLRINT(XB,YB,NB,IREAD,IDUF)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON/ModRes/ntr,xres,yres,eres,nrbin,ermin,ermax
      REAL XB(m_d),YB(m_d),DER2(m_d),xr(m_d),yr(m_d),er(m_d)
      real xres(m_d),yres(m_d),eres(m_d)
      LOGICAL  LSTART
      DATA     LSTART /.FALSE./
      if(IREAD.eq.0)then
       LSTART=.FALSE.
      else
       LSTART=.TRUE.
      endif
      SMALL=1.0E-20
      NFFT=m_d
      DO I=1,NB
       xr(I)=xres(I)
       yr(I)=yres(I)
       er(I)=eres(I)
      END DO
      CALL BINBLR(xr,yr,er,NB,XB,YB,NRBIN)
      XDMIN=XDAT(1)
      XDMAX=XDAT(NDAT)
      YMAX=0.0
      YSUM=0.0
      DO 10 I=1,NB
        IF (XB(I).LT.XDMIN .OR. XB(I).GT.XDMAX) GOTO 10
        IF (YB(I).GT.YMAX) YMAX=YB(I)
        YSUM=YSUM+YB(I)
  10  CONTINUE
      IF (YSUM.LT.SMALL) THEN
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       write(53,1000)
1000   format(' Ysum is too small')
       IDUF=1
       close(unit=53)
       RETURN
      ENDIF
      CALL XGINIT(XB,YB,NB,YMAX,XBMIN,XBMAX,LSTART)
      CALL SPLINE(XB,YB,NB,0.0,0.0,DER2)
      TWOPIN=2.0*3.141592654/FLOAT(NFFT)
      CALL VRFILL(FRES,0.0,NFFT)
      XX=0.0
      DXJ=XJ(2)-XJ(1)
      CALL SPLINT(XB,YB,DER2,NB,XX,FRES(1))
      SUM=FRES(1)
      do I=1,NFFT/2
        XX=XX+DXJ
        IF (XX.LT.XBMAX) CALL SPLINT(XB,YB,DER2,NB,XX,FRES(I+1))
        IF (-XX.GT.XBMIN) CALL SPLINT(XB,YB,DER2,NB,-XX,FRES(NFFT+1-I))
        SUM=SUM+FRES(I+1)+FRES(NFFT+1-I)
        TWOPIK(I)=TWOPIN*FLOAT(I-1)
      end do
      TWOPIK(NFFT/2+1)=TWOPIN*FLOAT(NFFT/2)
      BNORM=1.0/(SUM*FLOAT(NFFT))
      do I=1,NFFT
        FRES(I)=BNORM*FRES(I)
      end do
      CALL FOUR2(FRES,NFFT,1,1,0)
      do I=3,NFFT,4
        FRES(I)=-FRES(I)
        FRES(I+1)=-FRES(I+1)
      end do
      IF (.NOT. LSTART) THEN
        CALL VCOPY(FRES,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
      ENDIF
      LSTART=.TRUE.
      END
C     -----------------------------------------
      SUBROUTINE BINBLR(WX,WY,WE,NB,XB,YB,NBIN)
      INCLUDE 'mod_files.f90'
      REAL WX(*),WY(*),WE(*),XB(*),YB(*)
      N=0
      SMALL=1.0E-20
      BNORM=1.0/FLOAT(NBIN)
      do I=1,NB,NBIN
        N=N+1
        XXD=0.0
        DD=0.0
        K=0
        do J=0,NBIN-1
         IJ=I+J
          IF (IJ.LE.NB) THEN
            XXD=XXD+WX(IJ)
            IF (WE(IJ).GT.SMALL) THEN
              K=K+1
              DD=DD+WY(IJ)
            ENDIF
          ENDIF
        end do
        XB(N)=BNORM*XXD
        YB(N)=0.0
        IF (K.GT.0) YB(N)=BNORM*DD
      end do
      NB=N
      close(unit=53)
      END
C     ----------------------------------------------------------------
      SUBROUTINE XGINIT(XB,YB,NB,YMAX,XMIN,XMAX,LST)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      REAL    XB(*),YB(*)
      logical LST
      XDMIN=XDAT(1)
      XDMAX=XDAT(NDAT)
      Y0=YMAX/10.0
      do I=1,NB
       IF (YB(I).GE.Y0 .AND. XB(I).GT.XDMIN) GOTO 1
      end do
   1  XMIN=XB(I)
      do I=NB,1,-1
       IF (YB(I).GE.Y0 .AND. XB(I).LT.XDMAX) GOTO 2
      end do
   2  XMAX=XB(I)
      BWIDTH=XMAX-XMIN
      DXJ=BWIDTH/20.0
      AXMAX=ABS(XDAT(1))
      IF (ABS(XDAT(NDAT)).GT.AXMAX) AXMAX=ABS(XDAT(NDAT))
      XNDMAX=500.0
      IF (NDAT.GT.INT(XNDMAX)) XNDMAX=FLOAT(NDAT)
      DXDAT=2.0*AXMAX/XNDMAX
      IF (DXDAT.GT.DXJ) DXJ=DXDAT
      XNGD=(2.0*AXMAX)/DXJ
      NGD=NINT(LOG(XNGD-1.0)/LOG(2.0))+1
      NGD=2**NGD
	  IF (NGD.GT.m_d) then
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       WRITE(53,*)' ERROR in XGINIT : too many points'
       close(unit=53)
       RETURN
      ENDIF
      NFFT=NGD
      XJ(1)=-DXJ*FLOAT(NFFT/2)
      do I=2,NFFT
       XJ(I)=XJ(I-1)+DXJ
      end do
      XMIN=XMIN-5.0*BWIDTH
      XMAX=XMAX+5.0*BWIDTH
      IF (XMIN.LT.XB(1)) XMIN=XB(1)
      IF (XMAX.GT.XB(NB)) XMAX=XB(NB)
      IF (XMIN.LT.XDMIN) XMIN=XDMIN
      IF (XMAX.GT.XDMAX) XMAX=XDMAX
      if(LST)then
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1  access='append')
       WRITE(53,110) XMIN,XMAX
110    FORMAT(' Resolution Range: ',F7.1,'  to ',F7.1,'  ueV')
       close(unit=53)
      endif
      do I=1,NB
       IF (XB(I).GE.XMIN) GOTO 4
      end do
   4  IMIN=I
      do I=NB,1,-1
       IF (XB(I).LE.XMAX) GOTO 5
      end do
   5  IMAX=I
      B1=0.0
      B2=0.0
      do I=1,5
        B1=B1+YB(IMIN+I-1)
        B2=B2+YB(IMAX-I+1)
      end do
      B1=B1/5.0
      B2=B2/5.0
      DB=(B2-B1)/FLOAT(MAX(IMAX-IMIN-4,1))
      B=B1
      do I=IMIN,IMAX
        YB(I)=YB(I)-B
        B=B+DB
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
      REAL    V(*)
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
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      REAL      GRAD(*),HESS(*),DPAR(*),COVAR(*),FITP(*)
      INTEGER   INDX(*)
      if(o_w1.eq.1)then
       if (NFEW.GE.1) then
          FITP(5)=0.1
          FITP(6)=QW1(ISPEC)/WSCL
          if (NFEW.EQ.1) then
           CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
           RETURN
          endif
       endif
      endif
      J=4+2*NFEW
      DXLOG=0.85
      NSRCH=NINT(LOG(5.0*GSCL/WSCL)/LOG(DXLOG))
      CMIN=1.0E20
      FITP(J-1)=0.1
      FITP(J)=1.0
      do I=1,NSRCH
        CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
        CNORM=CCHI(FITP)
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
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      REAL            GRAD(*),HESS(NP,*),COVAR(NP,*)
      INTEGER         INDX(*)
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
      do I=1,NFEW
        J=3+I+I
        AJ=FITP(J)*ASCL
        CALL VMLTRC(EXPF(1,I),FR2PIK,NFT2,WORK)
        CALL VCOPY(WORK,FWRK,NFFT+2)
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,DDDPAR(1,J))
        CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        CALL VCOPY(FWRK,WORK,NFFT+2)
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
        CALL FOUR2(FWRK,NFFT,1,-1,-1)
        CALL DEGRID(FWRK,WORK(1,2))
        CALL VRDOTR(RESID,WORK(1,2),NDAT,SM)
        HESS(4,J+1)=-AJ*SM
        HESS(J+1,4)=HESS(4,J+1)
        CALL VMLTIC(WORK,NFT2,WORK)
        CALL FOUR2(WORK,NFFT,1,-1,-1)
        CALL DEGRID(WORK,FWRK)
        CALL VRDOTR(RESID,FWRK,NDAT,SM)
        HESS(J+1,J+1)=-AJ*SM
      end do
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC(1,2))
      CALL HESS1(HESS,NP,SCLVEC(1,2),STEPSZ,NFEW)
      if(o_w1.eq.1)then
        if (NP.GE.6) then
          DIF=WSCL*FITP(6)-QW1(ISPEC)
          SIG2=2.0/SIGQW1(ISPEC)**2
          GRAD(6)=GRAD(6)+SIG2*DIF*WSCL
          HESS(6,6)=HESS(6,6)+SIG2*WSCL**2
        endif
      endif
      CALL INVERT(HESS,COVAR,NP,INDX,DETLOG)
      END
C***<see the fit>*******************************************************
C
      SUBROUTINE SEEFIT(SIGPAR,CNORM,PRMSV,SIGSV)
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
      WRITE(53,150) FITP(4)*GSCL*1000.0
     1 ,SIGPAR(4)*GSCL*ERRSCL*1000.0
 150  FORMAT(' Zero offset       = ',F12.2,'  +- ',F10.2,'  ueV')
      WRITE(53,*)' Elastic line'
      WRITE(53,160) FITP(3)*ASCL,SIGPAR(3)*ASCL*ERRSCL
 160  FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
      PRMSV(1)=FITP(3)*ASCL
      SIGSV(1)=SIGPAR(3)*ASCL*ERRSCL
      do I=1,NFEW
        J=4+I+I
        WRITE(53,161) I
 161    FORMAT(' Quasi-elastic line ',I2)
        WRITE(53,170) 2000.0*FITP(J)*WSCL,
     1   2000.0*SIGPAR(J)*WSCL*ERRSCL
 170    FORMAT(5X,' FWHM        =   ',F13.2,'  +- ',F11.2,' ueV')
        WRITE(53,180) FITP(J-1)*ASCL,SIGPAR(J-1)*ASCL*ERRSCL
 180    FORMAT(5X,' Amplitude  =   ',1PE13.4,'  +- ',E11.3)
        PRMSV(I+I)=FITP(J-1)*ASCL
        SIGSV(I+I)=SIGPAR(J-1)*ASCL*ERRSCL
        PRMSV(I+I+1)=2.0*FITP(J)*WSCL
        SIGSV(I+I+1)=2.0*SIGPAR(J)*WSCL*ERRSCL
      end do
      close(unit=53)
      END
C
C***<numerical recipes routines>****************************************
C
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=10000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      do I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      end do
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      do K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
      end do
      RETURN
      END
C
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
