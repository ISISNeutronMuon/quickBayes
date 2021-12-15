C  PROGRAM RESNORM
C-----------------------------------------------------------------------
C  This is a simple resolution normalisation program, which was modified
C  from QUASI_LINES_2D.FOR . It reads in a 2-D resolution file, which 
C  should be of the same intermediate GENIE binary format as that of the 
C  2-D quasi-elastic data to be subsequently analysed. The output is in
C  file FOR007.DAT .
C-----------------------------------------------------------------------
C  Written by: D.S. Sivia, Rutherford Appleton Lab., Oxfordshire, England.
C-----------------------------------------------------------------------
C  Initial release for trials by CJC, MAA, WSH.       (DSS: 10-MAR-1992)
C  Stretch-factor range increased to +/- 10%.         (DSS: 27-MAY-1992)
C  WBF array-size increased in subroutine INBUFF.     (DSS: 28-MAY-1992)
C  Also put in READONLY for the resolution file.      (DSS: 28-MAY-1992)
C  Fixed bug in CLRBUF; problem for lots of detectors.(DSS: 29-MAY-1992)
C  Hard-wired meV input & ueV output, increased search 
C  to +/- 15% and added silent-running option.        (DSS:  2-OCT-1992)
C  Took out small and negative number safety-valve.   (DSS: 17-DEC-1992)
C  Fixed bug in DATIN1 (no more large error-bars!).   (DSS: 27-SEP-1993)
C-----------------------------------------------------------------------
c	nd,xout,yout,eout,yfit=resnorm.resnorm(numb,Xv,Yv,Ev,reals,Xdat,Xb,Yb,wrk,lwrk)
c	numb = [ngrp, nsp, ntc, Ndat, nbin, Imin, Imax, Nb]
c	reals = [efix, theta[0], rscl, bnorm]
c
      SUBROUTINE RESNORM(numb,x_in,y_in,e_in,reals,
     1 XD_in,XB_in,YB_in,vfile,rfile,l_fn,
     2 nd_out,xout,yout,eout,yfit,pfit)
      INCLUDE 'mod_data.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      common/mess/ness
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /ModPars/ NBIN,IMIN,IMAX,RSCL,BNORM
      real x_in(m_d), y_in(m_d), e_in(m_d)
cf2py intent(in) :: x_in, y_in, e_in                    !sample data
      real reals(4)
cf2py intent(in) :: reals                               !real parameters
      real XD_in(m_d), XB_in(m_d), YB_in(m_d)
cf2py intent(in) :: XD_in, XB_in, YB_in                 !sample xrange, res data (blur)
      integer numb(8)
cf2py intent(in) :: numb                                !integer parameters
      integer l_fn
cf2py intent(in) :: l_fn                                 !length of filenames
      character*140 vfile, rfile
cf2py intent(in) :: vfile, rfile                         !van & res filenames
      integer nd_out
cf2py intent(out) :: nd_out                              !number of output points
      real xout(m_d), yout(m_d), eout(m_d), yfit(m_d), pfit(4)
cf2py intent(out) :: xout, yout, eout, yfit, pfit              !data values
      real*4 XBLR(m_d),YBLR(m_d)
      real*4 GRAD(m_p),HESS(m_p,m_p),DPAR(m_p),COVAR(m_p,m_p)
      real*4 SIGPAR(m_p),FITPSV(m_p)
      real*4 PRBSV(4,m_sp),PRMSV(7,4,m_sp),SIGSV(7,4,m_sp)
      real*4 CHIX(31),BPARSV(m_p),CCHIX(31)
      real*4 DTNORM(m_sp),XSCALE(m_sp)
      integer INDX(m_p)
      character*80 text,TITLE
      DATA NCHIX /31/
      prog='r'
      o_bgd=1
      o_w1=0
      CALL init_paras
      do n=1,m_p
       GRAD(n)=0.0
       DPAR(n)=0.0
       SIGPAR(n)=0.0
       FITPSV(n)=0.0
       BPARSV(n)=0.0
       do m=1,m_p
        HESS(n,m)=0.0
        COVAR(n,m)=0.0
       end do
      end do
      do n=1,m_sp
       do n1=1,4
        PRBSV(n1,n)=0.0
        do n2=1,7
         PRMSV(n2,n1,n)=0.0
         SIGSV(n2,n1,n)=0.0
        end do
       end do
      end do
      do n=1,31
       CHIX(n)=0.0
       CCHIX(n)=0.0
      end do
c
c     numb = [ngrp, nsp, ntc, Ndat, nbin, Imin, Imax, Nb]
      NSPEC=numb(1)                                      !no. of groups
      ISP=numb(2)                                   !group number
      ntc=numb(3)                                     !no. of points
      NDAT=numb(4)
      NBIN=numb(5)
      IMIN=numb(6)
      IMAX=numb(7)
      NB=numb(8)
c     reals = [efix, theta[isp], rscl, bnorm]
      efix=reals(1)
      theta(ISP)=reals(2)
      RSCL=reals(3)
      BNORM=reals(4)
      do n=1,m_d
       xin(n)=x_in(n)
       yin(n)=y_in(n)
       ein(n)=e_in(n)
       XDAT(n)=XD_in(n)
       XBLR(n)=XB_in(n)
       YBLR(n)=YB_in(n)
      end do
      do n=1,m_sp                          !ResNorm output not used
       DTNORM(n)=1.0                       !DTNORM, NSPEC
       XSCALE(n)=1.0                       !XSCALE, NSPEC
      end do
	    lptfile = ''
      l_lpt=l_fn+11
      lptfile(1:l_lpt)=vfile(1:l_fn)//'_resnrm.lpt'
      if(ISP.eq.1)then                  !print info	
       call open_f(53,lptfile)
       write(53,1107)vfile
1107   format(' Vanadium file : ',a140)
       write(53,1108)rfile
1108   format(' Resolution file : ',a140)
       close(unit=53)
      endif
      NDAT=ntc
      CALL BLRINT(XBLR,YBLR,NB,1.0,.false.)
      CALL DPINIT
      CALL GDINIT
      DCHIX=0.010
      CHIX0=1.0-DCHIX*FLOAT((NCHIX-1)/2)
c spectrum starts
      NDAT=ntc
      CALL DATIN(ISP,DTNORM,IDUF)
      if(ISP.eq.1)then                  !print info	
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       write(53,1110)xin(1),xin(ntc)
1110   format(' Energy range: ',f8.3,' to ',f8.3,' meV')
       write(53,1002)xblr(1),xblr(NB+1)
1002   format(' Energy range after BLRint: ',f8.3,' to ',f8.3,' meV')
       close(unit=53)
      endif
      CMINSV=1.0E+20
      DO IX=1,NCHIX
        CHIX(IX)=CHIX0+FLOAT(IX-1)*DCHIX
        CALL BLRINT(XBLR,YBLR,NB,CHIX(IX),.false.)
        CALL DPINIT
        CALL PRINIT(FITP,3,NFEW,IX)
        CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
        NPARMS=4+2*NFEW
        CHIOLD=CCHI(FITP)
        CALL VCOPY(FITP,FITPSV,NPARMS)
        STEPSZ=0.3
        IF (NFEW.GT.1) STEPSZ=STEPSZ/10.0
        IAGAIN=0
        CDIFMN=0.003
        ness=0
        DO I=1,200
          CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,
     1    STEPSZ,.FALSE.)
          CALL NEWEST(COVAR,GRAD,NPARMS,NFEW,DPAR,FITP)
          CNORM=CCHI(FITP)
          IF (CNORM.LE.CHIOLD) THEN
            CHIDIF=(CHIOLD-CNORM)/CNORM
            IF (ABS(CHIDIF).LE.CDIFMN) THEN
              IF (IAGAIN.EQ.0) THEN
                CDIFMN=0.00005
                STEPSZ=0.15
                IAGAIN=1
              ELSE
                GOTO 3
              ENDIF
            ENDIF
            CHIOLD=CNORM
            CALL VCOPY(FITP,FITPSV,NPARMS)
          ELSE
            CALL VCOPY(FITPSV,FITP,NPARMS)
            STEPSZ=STEPSZ*0.6
          ENDIF
        END DO
   3    CCHIX(IX)=CNORM
        IF (CNORM.LT.CMINSV) THEN
          IXMIN=IX
          CMINSV=CNORM
          XSCLSV=CHIX(IX)
          FLUXSV=FITP(3)*ASCL
          CALL VCOPY(FITP,BPARSV,NPARMS)
        ENDIF
      END DO
      CALL VCOPY(BPARSV,FITP,NPARMS)
      CALL QUADRT(CCHIX,IXMIN,NCHIX,XSCLSV,CHIX)
      CALL BLRINT(XBLR,YBLR,NB,XSCLSV,.false.)
      CALL DPINIT
      CALL PRINIT(FITP,3,NFEW,2)
      do I=1,5
        CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,
     1  0.01,.TRUE.)
        CALL NEWEST(COVAR,GRAD,NPARMS,NFEW,DPAR,FITP)
      end do
      CALL ERRBAR(COVAR,NPARMS,SIGPAR)
      CNORM=CCHI(FITP)
      pfit(1)=FITP(3)*ASCL
      pfit(2)=XSCLSV
      pfit(3)=CNORM
      pfit(4)=FITP(4)*GSCL
      OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
      write(53,110) pfit(1),pfit(2),pfit(3),pfit(4)
 110  FORMAT(G15.5,F8.3,5X,2F10.3)
      close(unit=53)

      nd_out=NDAT
      do n=1,nd_out
       xout(n)=XDAT(n)
       yout(n)=DAT(n)
       if(SIG(n).gt.1.0e-20)then
        eout(n)=SQRT(2.0/SIG(n))
       else
        eout(n)=0.0
       endif
       yfit(n)=FIT(n)
      end do
c spectrum ends
      CLOSE(unit=53)
      RETURN
      END
