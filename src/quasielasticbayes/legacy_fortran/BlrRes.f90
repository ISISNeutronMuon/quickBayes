C***<set up blur function>**********************************************
      SUBROUTINE BLRINT(XB,YB,NB,XSCALE,LSTART)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      REAL XB(m_d),YB(m_d)
      LOGICAL LSTART
      NFFT=m_d
      IF (.NOT .LSTART) THEN
       YMAX=0.0
       do I=1,NFFT
        IF (YB(I).GT.YMAX) YMAX=YB(I)
       end do
      ELSE
       OPEN(UNIT=53,FILE=lptfile,STATUS='old',FORM='formatted',
     1 access='append')
       WRITE(53,*)' Resolution function X-axis multiplied by: ',XSCALE
       close(unit=53)
      ENDIF
      CALL XGINIT(XB,YB,NB,YMAX,XSCALE)
      DXJ=XJ(2)-XJ(1)
      DXB=XB(2)-XB(1)
      do I=1,NB
       IF (XB(I).GE.0.0) GOTO 4
      end do
   4  IR0=I
      XB0=XB(IR0)
      DXB=DXB*XSCALE
      XBMIN=XB(1)*XSCALE+DXB
      XBMAX=XB(NB)*XSCALE-DXB
      TWOPIN=2.0*3.141592654/FLOAT(NFFT)
      CALL VRFILL(FRES,0.0,NFFT)
      X=0.0
      CALL RINTRP(X,YB,IR0,XB0,DXB,FRES(1))
      SM=FRES(1)
      do I=1,NFFT/2
       X=X+DXJ
       IF (X.LT.XBMAX) CALL RINTRP(X,YB,IR0,XB0,DXB,FRES(I+1))
       IF (-X.GT.XBMIN) CALL RINTRP(-X,YB,IR0,XB0,DXB,FRES(NFFT+1-I))
       SM=SM+FRES(I+1)+FRES(NFFT+1-I)
       TWOPIK(I)=TWOPIN*FLOAT(I-1)
      end do
      TWOPIK(NFFT/2+1)=TWOPIN*FLOAT(NFFT/2)
      BNORM=1.0/(SM*FLOAT(NFFT))
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
      END
C     ---------------------------------------------------------
      SUBROUTINE XGINIT(XB,YB,NB,YMAX,XSCALE)
      INCLUDE 'res_par.f90'
      INCLUDE 'mod_files.f90'
      COMMON /FFTCOM/ FRES(m_d2),FWRK(m_d2),XJ(m_d),TWOPIK(m_d1),NFFT
      COMMON /DATCOM/ XDAT(m_d),DAT(m_d),SIG(m_d),NDAT
      REAL XB(*),YB(*)
      Y0=YMAX/10.0
      do I=1,NB
       IF (YB(I).GE.Y0) GOTO 1
      end do
   1  IMIN=I
      do I=NB,1,-1
       IF (YB(I).GE.Y0) GOTO 2
      end do
   2  IMAX=I
      BWIDTH=(XB(IMAX)-XB(IMIN))*XSCALE
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
      END
C     --------------------------------------
      SUBROUTINE RINTRP(X,YB,IR0,XB0,DX,BLR)
      REAL YB(*)
      XOFSET=(X-XB0)/DX
      IXOFST=INT(XOFSET)
      DXOFST=ABS(XOFSET-FLOAT(IXOFST))
      I=IXOFST+IR0
      IF (XOFSET.GE.0.0) THEN
       BLR=YB(I)+DXOFST*(YB(I+1)-YB(I))
      ELSE
       BLR=YB(I)+DXOFST*(YB(I-1)-YB(I))
      ENDIF
      END
