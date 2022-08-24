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
C     ----------------------------------------
      SUBROUTINE QUADRT(CHI,IXMIN,NX,XSCLSV,X)
      real*4 CHI(*),X(*)
      IF (IXMIN.LE.1 .OR. IXMIN.GE.NX) RETURN
      CM=CHI(IXMIN-1)
      C0=CHI(IXMIN)
      CP=CHI(IXMIN+1)
      X0=0.5*(CM-CP)/(CM-C0-C0+CP)
      IF (ABS(X0).LT.0.5) XSCLSV=X0*(X(2)-X(1))+X(IXMIN)
      END
c
C***<calculate data & chi-squared>**************************************
C
      FUNCTION CCHI(V)
      INCLUDE 'res_par.f90'
      common/mess/ness
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      COMMON /WRKCOM/ WORK(m_d2,2)
      real*4 V(*)
      CHI=0.0
      B1=BSCL*V(1)
      B2=BSCL*V(2)
      A0=ASCL*V(3)
      DELTAX=V(4)
      NFT2=NFFT/2+1
      CALL CXSHFT(FRES,DELTAX,TWOPIK,FR2PIK,FR2PIK(1,2),NFT2)
      CALL VRFILL(WORK,A0,NFT2)
      XNSCL=-LOG(1.0E-7)/(TWOPIK(2)-TWOPIK(1))
      DO J=1,NFEW
        AJ=ASCL*V(3+J+J)
        SIGJ=WSCL*V(4+J+J)/GSCL
        CALL VRFILL(EXPF(1,J),0.0,NFT2)
        NXFT2=1+NINT(XNSCL/(ABS(SIGJ)+1.0E-10))
        IF (NXFT2.GT.NFT2) NXFT2=NFT2
        DO I=1,NXFT2
          EXPIJ=EXP(-TWOPIK(I)*SIGJ)
          EXPF(I,J)=EXPIJ
          WORK(I,1)=WORK(I,1)+AJ*EXPIJ
        END DO
      END DO
      CALL VMLTRC(WORK,FR2PIK,NFT2,FWRK)
      CALL FOUR2(FWRK,NFFT,1,-1,-1)
      CALL DEGRID(FWRK,FIT)
      X1=XDAT(1)
      BNORM=(B2-B1)/(XDAT(NDAT)-X1)
      DO I=1,NDAT
        FIT(I)=FIT(I)+B1+BNORM*(XDAT(I)-X1)
        DIF=FIT(I)-DAT(I)
        RESID(I)=DIF*SIG(I)
        CHI=CHI+DIF*RESID(I)
      END DO
      CCHI=CHI/(2.0*FLOAT(NDAT))
      END
C
C**<search for one more & refine amplitudes>****************************
C
      SUBROUTINE REFINE(GRAD,HESS,NP,DETLOG,INDX,COVAR,STEPSZ,lpt)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      COMMON /DINTRP/ IPDAT(m_d),XPDAT(m_d)
      COMMON /FITCOM/ FIT(m_d),RESID(m_d),NFEW,FITP(m_p),EXPF(m_d1,6)
      COMMON /SCLCOM/ BSCL,ASCL,WSCL,SCLVEC(m_p,2),GSCL
      COMMON /WRKCOM/ WORK(m_d2,2)
      COMMON /GRDCOM/ DDDPAR(m_d,m_p),FR2PIK(m_d2,2)
      real*4 GRAD(*),HESS(NP,*),COVAR(NP,*)
      integer INDX(*)
      logical lpt
      NFT2=NFFT/2+1
      CNORM=CCHI(FITP)
      if(lpt)then
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       WRITE(53,100) CNORM
 100   FORMAT(' Refine: Normalised chi-squared = ',1PE13.5)
       close(unit=53)
      endif
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
      DO I=1,NFEW
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
      END DO
      CALL GRADPR(GRAD,RESID,NDAT,NP,SCLVEC(1,2))
      CALL HESS1(HESS,NP,SCLVEC(1,2),STEPSZ,NFEW)
      CALL INVERT(HESS,COVAR,NP,INDX,DETLOG)
      END
