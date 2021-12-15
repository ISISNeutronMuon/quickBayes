c      PROGRAM QUASI_STEXP_2DPLOT0
C     ---------------------------
C
C-----------------------------------------------------------------------
C  This program tries to estimate the parameters of an elstic peak and 
C  a single "stretched exponential" pertaining to some relevant data.
C  It was modified from QUASI_LINES_2D, a 2-D Bayesian Quasi-elastic 
C  line-fitting program. The resolution function is given on a constant 
C  X-binning grid, and is assumed to be invariant. The background is 
C  assumed to be linear and the data should be of an intermediate Genie
C  binary format, such as that output by batch-mode ICON.
C-----------------------------------------------------------------------
C  Written by: D.S. Sivia, Rutherford Appleton Lab., Oxfordshire, England.
C-----------------------------------------------------------------------
C Initial release for trials by CJC, MAA, WSH.        (DSS: 27-FEB-1992)
C Changed file access from sequential to direct.      (DSS:  5-MAR-1992)
C Changed STEPSZ reduction from 0.85 to 0.6 .         (DSS:  7-MAR-1992)
C Can read in a detector normalisation file.          (DSS: 12-MAR-1992)
C Readonly for resol. file; silent-running modified.  (DSS: 14-MAY-1992)
C WBF array-size increased in subroutine INBUFF.      (DSS: 28-MAY-1992)
C Fixed bug in CLRBUF; problem for lots of detectors. (DSS: 29-MAY-1992)
C Hard-wired meV as input and ueV as output.          (DSS:  2-OCT-1992)
C Put in STEPSZ safety-valve in case Chisq stuck!     (DSS: 27-OCT-1992)
C Took out small and negative number safety-valve.    (DSS: 17-DEC-1992)
C Fixed bug in DATIN1 (no more large error-bars!).    (DSS: 27-SEP-1993)
C Modification of QUASI_LINES_2D.FOR .                (DSS: 26-APR-1994)
C Do "stretched-exp" quasi_lines properly!            (DSS:  7-JUN-1999)
C Plot the 2d posterior pdf for Beta and Sigma.       (DSS: 10-JUN-1999)
C-----------------------------------------------------------------------
C
      SUBROUTINE Quest(numb,x_in,y_in,e_in,reals,opft,
     1 XD_in,XB_in,YB_in,sfile,rfile,l_fn,
     2 xs_out,ys_out,xb_out,yb_out,zp_out)
      INCLUDE 'mod_data.f90'
      INCLUDE 'mod_files.f90'
      INCLUDE 'options.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /QW1COM/ QW1(m_sp),SIGQW1(m_sp),ISPEC
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /STEXP/  BETEXP,SIGMA,ADETLG
      COMMON/ModPars/NBIN,IMIN,IMAX,RSCL,BNORM
      integer numb(11)
cf2py intent(in) :: numb                                !integer parameters
      real x_in(m_d), y_in(m_d), e_in(m_d)
cf2py intent(in) :: x_in, y_in, e_in                    !sample data
      real reals(4)
cf2py intent(in) :: reals                               !real parameters
      integer opft(4)
cf2py intent(in) :: opft                                 !options parameters
      real XD_in(m_d), XB_in(m_d), YB_in(m_d)
cf2py intent(in) :: XD_in, XB_in, YB_in                 !sample xrange, res data (blur)
      character*140 sfile, rfile
cf2py intent(in) :: sfile, rfile                         !sample & res filenames
      integer l_fn
cf2py intent(in) :: l_fn                                 !length of filenames
      real xs_out(201), ys_out(201), xb_out(201), yb_out(201)
cf2py intent(out) :: xs_out, ys_out, xb_out, yb_out
      real zp_out(40500)
cf2py intent(out) :: zp_out
      integer l_title,l_user
      character*80 user,title
      real XBLR(m_d),YBLR(m_d)
      REAL GRAD(m_p),HESS(m_p,m_p),DPAR(m_p),COVAR(m_p,m_p)
      REAL SIGPAR(m_p),FITPSV(m_p)
      REAL DTNORM(m_sp),XSCALE(m_sp)
      INTEGER INDX(m_p),NBMIN(2)
      REAL XBETA(201),YSIGMA(201,2),PRBETA(201),PRSIGM(201),xsig(201)
      REAL      ZP(40500),ZLP(40500)
      LOGICAL   LBATCH
      DATA      LBATCH,NBET,NSIG /.FALSE.,201,201/
      prog='q'
      CALL init_paras
C
c     numb = [ngrp, nsp, ntc, Ndat, nbin, Imin, Imax, Nb, nrbin, Nbet,Nsig]
      NSPEC=numb(1)                                      !no. of groups
      ISP=numb(2)                                   !group number
      ISPEC=ISP
      ntc=numb(3)                                     !no. of points
      NDAT=numb(4)
      NBIN=numb(5)
      IMIN=numb(6)
      IMAX=numb(7)
      NB=numb(8)
      nrbin=numb(9)
      NBET=numb(10)
      NSIG=numb(11)
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
      o_el=opft(1)
      o_bgd=opft(2)
      o_w1=opft(3)
      do n=1,m_sp
       DTNORM(n)=1.0                       !'dtnorm', DTNORM, NSPEC
       XSCALE(n)=1.0                       !'xscale', XSCALE, NSPEC
       yprob=0.0
      end do
      lptfile = ''
      fileout1 = ''
      fileout2 = ''
      fileout3 = ''
      l_lpt=l_fn+8
      lptfile(1:l_lpt)=sfile(1:l_fn)//'_Qst.lpt'
      l_file=l_fn+8
      fileout1(1:l_lpt)=sfile(1:l_fn)//'_Qsb.ql1'
      fileout2(1:l_lpt)=sfile(1:l_fn)//'_Qss.ql2'
      l_title=9
      title='<unknown>'
      l_user=9
      user='<unknown>'
      o_w1=0
c
      if(ISP.eq.1)then                      !print info	
       call open_f(53,lptfile)
       write(53,1107)sfile
1107   format(' Sample file : ',a140)
       write(53,1108)rfile
1108   format(' Resolution file : ',a140)
       write(53,1110)xin(imin),xin(imax)
1110   format(' Energy range: ',f8.3,' to ',f8.3,' meV')
       if (o_el.eq.0) write(53,1111)
       if (o_el.eq.1) write(53,1112)
       if (o_bgd.eq.2) write(53,1113)
       if (o_bgd.eq.1) write(53,1114)
       if (o_bgd.eq.0) write(53,1115)
1111   format(' Elastic option : NO peak')
1112   format(' Elastic option : WITH peak')
1113   format(' Background option : sloping')
1114   format(' Background option : flat')
1115   format(' Background option : zero')
       if(o_w1.eq.1)then
        WRITE(53,*)'Width option : fixed from file ',difile
       else
        WRITE(53,*)'Width option : free '
       endif
       close(unit=53)
       call open_f(1,fileout1)
       call open_f(2,fileout2)
       do n=1,2
        write(n,1107)sfile
        write(n,1121)title(1:l_title)
1121    format(' Title : ',a)
        write(n,1122)user(1:l_user)
1122    format(' User : ',a)
        if(n.eq.1)WRITE(n,120)NSPEC,NBET
        if(n.eq.2)WRITE(n,120)NSPEC,NSIG
 120    FORMAT(2X,3I10)
        WRITE(n,*)' -------------------------------------------------'
        close(unit=n)
       end do
      endif                             !end print
c
      CALL BLRINT(XBLR,YBLR,NB,1.0,.false.)
      CALL DPINIT
      CALL GDINIT
      CALL BTINIT(XBETA,NBET,0.1,1.4,NBMIN)
        BETEXP=1.0
        ISPEC=ISP
        CALL DATIN(ISP,DTNORM,IDUF)
        CALL BLRINT(XBLR,YBLR,NB,XSCALE(ISP),.false.)
        CALL DPINIT
        CALL PRINIT(FITP,3,NFEW,1)
        CALL REFINA(GRAD,HESS,DPAR,3+NFEW,DETLOG,INDX,COVAR)
   1    NPARMS=4+2*NFEW
        CHIOLD=CBCHI(FITP)
        CALL VCOPY(FITP,FITPSV,NPARMS)
        STEPSZ=0.3
        IF (NFEW.GE.1) STEPSZ=STEPSZ/10.0
        IAGAIN=0
        CDIFMN=0.003
        do I=1,200
          CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,STEPSZ)
          CALL NEWEST(COVAR,GRAD,NPARMS,NFEW,DPAR,FITP)
          CNORM=CBCHI(FITP)
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
            IF (STEPSZ.LT.1.0E-10) GOTO 3
          ENDIF
        end do
   3    IF (NFEW.EQ.0) THEN
          NFEW=NFEW+1
          CALL SEARCH(GRAD,HESS,DPAR,NFEW,INDX,COVAR,FITP)
          GOTO 1
        ENDIF
        CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.25)
        CNORM=CBCHI(FITP)
        CALL VCOPY(FITP,FITPSV,NPARMS)
        CALL REFINE(GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.7)
        CALL ERRBAR(COVAR,NPARMS,SIGPAR)
        SIG1SV=SIGPAR(4)
        CALL SGINIT(YSIGMA,201,NSIG,FITP(NPARMS),SIGPAR(NPARMS),WSCL)
        K=1
        ZLPMAX=-1.0E10
        do J=1,NSIG
          SIGMA=YSIGMA(J,1)
          do I=1,NBET
            BETEXP=XBETA(I)
            FITP(NPARMS+1)=BETEXP
            FITP(NPARMS)=SIGMA
            FITP(NPARMS-1)=0.1
            FITP(4)=FITPSV(1)
            CALL VCOPY(FITPSV,FITP,3)
            CALL REFINA(GRAD,HESS,DPAR,3+NFEW,ADETLG,INDX,COVAR)
            CHISQ=FLOAT(NDAT)*CBCHI(FITP)
            ZLP(K)=-0.5*(LOG10(2.7182818)*CHISQ+ADETLG)
            ZLPMAX=MAX(ZLPMAX,ZLP(K))
            K=K+1
          end do
        end do
        CALL PLNORM(ZLPMAX,ZLP,ZP,NBET*NSIG)
        CALL PRMARG(ZP,NBET,NSIG,PRBETA,PRSIGM)
      do n=1,NSIG
       xsig(n)=YSIGMA(n,2)
      end do
      OPEN(UNIT=1,FILE=fileout1,STATUS='old',FORM='formatted',
     1 access='append')
      write(1,130)(Xbeta(i),PRbeta(i),i=1,NBET)
130   format(2f10.5)
      close(unit=1)
      OPEN(UNIT=2,FILE=fileout2,STATUS='old',FORM='formatted',
     1 access='append')
      write(2,130)(xsig(i),PRsigm(i),i=1,NSIG)
      close(unit=2)
      nz=NBET*NSIG
      do n=1,Nbet
       xb_out(n)=XBETA(n)
       yb_out(n)=PRBETA(n)
      end do
      do n=1,Nsig
       xs_out(n)=XSIG(n)
       ys_out(n)=PRSIGM(n)
      end do
      do n=1,nz
       zp_out(n)=ZP(n)
      end do
      RETURN
      END
