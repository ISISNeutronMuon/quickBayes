# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from quasielasticbayes.python.fortran_python import *
from quasielasticbayes.python.constants import *
from quasielasticbayes.python.four import *
from quasielasticbayes.python.bayes import *
from math import pi
import numpy as np
from scipy.interpolate import interp1d
"""
<numerical recipes routines>****************************************
"""
def SPLINE(X,Y,N,YP1,YPN,Y2):
      """
      X,Y are vals
      YP1 is the derivative at start
      YPN is the '' at end
      Y2 is the list of second derivatives
      """
      return interp1d(X.output(), Y.output(), kind='cubic')

      NMAX=10000
      U = vec(NMAX)
      if YP1 > .99E30:
        Y2.set(1,0.)
        U.set(1,0.)
      else:
        Y2.set(1,-0.5)
        U.set(1,(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1))
      
      for I in get_range(2,N-1):
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2.set(I,(SIG-1.)/P)
        U.set(I,(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P)

      QN, UN = 0,0
      if YPN <= .99E30:
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))

      Y2.set(N,(UN-QN*U(N-1))/(QN*Y2(N-1)+1.))
      for K in get_range(N-1,1,-1):
        Y2.set(K,Y2(K)*Y2(K+1)+U(K))


def SPLINT(X, func):
      return func(X)

      KLO=1
      KHI=N
      while KHI-KLO <= 1:
        K=(KHI+KLO)/2
        if XA(K) > X:
          KHI=K
        else:
          KLO=K

      H=XA(KHI)-XA(KLO)
      if H==0:
         raise ValueError('Bad XA input.')
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

def rm_BG(x_bin,y_bin,N_bin, YMAX,LST,COMS, store, lptfile):
      # DATA XDAT NDAT, FFT XJ, N
      XDMIN=COMS["DATA"].XDAT(1)
      XDMAX=COMS["DATA"].XDAT(COMS["DATA"].NDAT)
      Y0=YMAX/10.0
      # these two check seem to be finiding the first/last y value greater than some ref value
      # that is also in the x range
      def check_data_from_below(XB, YB, Y0, XDMIN, I):
           return YB(I)>=Y0 and XB(I)>XDMIN
      I = find_index((x_bin, y_bin, Y0, XDMIN),1, N_bin, check_data_from_below)

      XMIN=x_bin(I)
      def check_data_from_above(XB, YB, Y0, XDMAX, I):
          return YB(I)>=Y0 and XB(I)<XDMAX
      I = find_index((x_bin, y_bin, Y0, XDMAX),N_bin,1, check_data_from_above, step=-1)
      XMAX=x_bin(I)

      # get x range of interesting values
      data_window = XMAX-XMIN
      bin_width =data_window/20.0

      AXMAX=abs(COMS["DATA"].XDAT(1))
      # assume data is approximatly symmetrical
      if abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))>AXMAX:
         AXMAX=abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))
      XNDMAX=500.0
      if COMS["DATA"].NDAT > int(XNDMAX):
         XNDMAX=float(COMS["DATA"].NDAT)

      mean_bin_width=2.0*AXMAX/XNDMAX # average bin width
      if mean_bin_width>bin_width: # use the largest bin width
         bin_width=mean_bin_width

      XNGD=(2.0*AXMAX)/bin_width # average bin size
      NGD=NINT(np.log(XNGD-1.0)/np.log(2.0))+1
      NGD=pow(2,NGD) # number of data points
      
      if NGD>m_d:
       store.open(53,lptfile)
       store.write(53,' ERROR in XGINIT : too many points')
       store.close(unit=53)
       return
      COMS["FFT"].NFFT=NGD
      COMS["res_data"].N_FT = NGD//2

      # Set the x values for evenly spaced bins
      COMS["FFT"].XJ.set(1, - bin_width*float(COMS["FFT"].NFFT/2))
      for j in get_range(2,COMS["FFT"].NFFT):
          COMS["FFT"].XJ.set(j,COMS["FFT"].XJ(j-1)+bin_width)
      # get the first and last bin
      XMIN=XMIN-5.0*data_window
      XMAX=XMAX+5.0*data_window
      if XMIN<x_bin(1):
         XMIN=x_bin(1)
      if XMAX>x_bin(N_bin):
         XMAX=x_bin(N_bin)
      if XMIN<XDMIN:
         XMIN=XDMIN
      if XMAX>XDMAX:
         XMAX=XDMAX
      if LST:
       store.open(53,lptfile)

       store.write(53,f' Resolution Range: {XMIN} to {XMAX} ueV')
       store.close(unit=53)
      
      # get x range -> via indices
      def check_data_x_min(XB, XMIN, I):
          return  XB(I)>=XMIN
      I = find_index((x_bin, XMIN),1,N_bin, check_data_x_min)
      IMIN=I

      def check_data_x_max(XB, XMAX, I):
          return  XB(I)<=XMAX
      I = find_index((x_bin, XMAX),N_bin, 1,check_data_x_max,step=-1)
      IMAX=I

      BG1=0.0
      BG2=0.0
      # get mean value for 5 bins (that are in range) closest to min/max
      for I in get_range(1,5):
        BG1=BG1+y_bin(IMIN+I-1)
        BG2=BG2+y_bin(IMAX-I+1)
      BG1=BG1/5.0
      BG2=BG2/5.0

      # the gradiant of the BG
      DB=(BG2-BG1)/float(max(IMAX-IMIN-4,1)) # no idea where the 4 comes from
      BG=BG1
      # Remove BG value from data
      for I in get_range(IMIN,IMAX):
        y_bin.set(I, y_bin(I)-BG)
        BG=BG+DB
      print("test",XMIN, XMAX, y_bin.output()[0:4] )
      return XMIN, XMAX, y_bin



@deprecated
def XGINIT(XB,YB,NB,YMAX,LST, COMS, store, lptfile):
      return rm_BG(XB,YB,NB,YMAX,LST, COMS, store, lptfile)
      # DATA XDAT NDAT, FFT XJ, N
"""
***<set up blur function>**********************************************
"""

def rebin(x_in,y_in,e_in,N_new_bins,N_merged_bins):
      """
      """
      XB = vec(N_new_bins)
      YB = vec(N_new_bins)
      SMALL=1.0E-20
      BNORM=1.0/float(N_merged_bins)

      for I in get_range(1,N_new_bins,N_merged_bins): # new binning
        x_value=0.0
        y_value=0.0
        K=0
        for J in get_range(0,N_merged_bins-1): # loop over bins in new bin
         IJ=I+J
         if IJ<=N_new_bins:

            x_value += x_in(IJ)
            if e_in(IJ) > SMALL: # only include non-zero errors
              K=K+1
              y_value += y_in(IJ)

         XB.set(I,BNORM*x_value)
         YB.set(I,0.0)
         if K>0:
            YB.set(I,BNORM*y_value) # normalise data
      return XB, YB

@deprecated
def BINBLR(WX,WY,WE,NB,NBIN):
      """
      Original dat is W* and the output is *B
      It seems to just be a rebin alg
      """
      return rebin(WX,WY,WE,NB,NBIN)


def bin_resolution(N_bin,IREAD,IDUF,COMS,store,lptfile):
      # data XDAT, NDAT

      LSTART = True
      if IREAD==0:
         LSTART=False

      SMALL=1.0E-20
      COMS["FFT"].NFFT=m_d
      COMS["res_data"].N_FT = m_d //2
      # copy resolution data
      xr = vec(m_d)
      yr = vec(m_d)
      er = vec(m_d)
      xr.copy(COMS["Res"].xres.output_range(1,N_bin))
      yr.copy(COMS["Res"].yres.output_range(1,N_bin))
      er.copy(COMS["Res"].eres.output_range(1,N_bin))

      # rebin the resolution data as it is not on an evenly spaced grid
      x_bin, y_bin = rebin(xr,yr,er,N_bin,COMS["Res"].nrbin)

      # sample data x range  
      XDMIN=COMS["DATA"].XDAT(1)
      XDMAX=COMS["DATA"].XDAT(COMS["DATA"].NDAT)

      # get indicies for valid x range
      first_index = np.argmax(x_bin.output() >=XDMIN)+1
      last_index = np.argmin(x_bin.output()<=XDMAX)+1
      # get total and max Y binned values within valid x range
      YSUM = np.sum(y_bin.output_range(first_index, last_index))
      YMAX = np.max(y_bin.output_range(first_index, last_index))
      if YSUM<SMALL:
       store.open(53,lptfile)
       store.write(53,' Ysum is too small')
       IDUF=1
       store.close(unit=53)
       return

      XBMIN, XBMAX, y_bin = rm_BG(x_bin,y_bin,N_bin,YMAX,LSTART, COMS,store, lptfile) # subtracts BG off y_bin

      # populate FRES with spline of evenly spaced binned data -> data to FFT later
      DER2 = vec(m_d)
      data = vec(COMS["FFT"].NFFT)
      func=SPLINE(x_bin,y_bin,N_bin,0.0,0.0,DER2)
      TWOPIN=np.pi/float(COMS["res_data"].N_FT)

      # clean up old data
      COMS["FFT"].FRES.fill(0.0, COMS["FFT"].NFFT)
      COMS["res_data"].FTY.fill(0.0, COMS["res_data"].N_FT)

      XX=0.0
      bin_width=COMS["FFT"].XJ(2)-COMS["FFT"].XJ(1)

      # spline the data ontop the sample data bins
      data = vec(COMS["FFT"].NFFT)
      data.set(1,SPLINT(XX,func))
      for I in get_range(1,int(COMS["FFT"].NFFT/2)):#use symmetry to have the range we need to loop
        XX += bin_width
        if XX < XBMAX:
           data.set(I+1,SPLINT(XX,func))
        if -XX > XBMIN:

           data.set(COMS["FFT"].NFFT+1-I,SPLINT(-XX,func))
        # set phases
        COMS["FFT"].TWOPIK.set(I,TWOPIN*float(I-1)) # looks to be the phase
        COMS["res_data"].phases.set(I,TWOPIN*float(I-1))
      COMS["FFT"].TWOPIK.set(int(COMS["FFT"].NFFT/2)+1,TWOPIN*float(COMS["FFT"].NFFT/2))
      COMS["res_data"].phases.set(int(COMS["res_data"].N_FT)+1,TWOPIN*float(COMS["res_data"].N_FT))
    
      # average y value in each bin -> integral of all data is 1
      SUM = np.sum(data.output())
      BNORM=1./(SUM*float(COMS["FFT"].NFFT))
      data.copy(data.output()*BNORM) # scale the splined data
      out = FOUR2(data, COMS["FFT"].NFFT,1,1,0)
      COMS["FFT"].FRES.copy(flatten(out))
      COMS["res_data"].FTY.copy(out)

      # some rotations?
      for I in get_range(3,COMS["FFT"].NFFT,4):
        COMS["FFT"].FRES.set(I,-COMS["FFT"].FRES(I))
        COMS["FFT"].FRES.set(I+1, -COMS["FFT"].FRES(I+1))

      for k in get_range(2, COMS["res_data"].N_FT, 2):
          COMS["res_data"].FTY.set(k, -COMS["res_data"].FTY(k))

      # IFT of resolution
      if not LSTART:
        tmp = COMS["FFT"].FRES.output_range(end=COMS["FFT"].NFFT+2)
        COMS["FFT"].FWRK.copy(tmp)
        out = FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
        COMS["FFT"].FWRK.copy(flatten(out[0:m_d2]))


        COMS["res_data"].IFTY.copy(COMS["res_data"].FTY.output_range(end=COMS["res_data"].N_FT))
        out2 = FOUR2(COMS["res_data"].IFTY, COMS["res_data"].N_FT, 1, -1,-1)
        COMS["res_data"].IFTY.copy(out[0:m_d2])

      LSTART= True
      return x_bin, y_bin
      #return x_bin, y_bin

@deprecated
def BLRINT(NB,IREAD,IDUF,COMS,store,lptfile):
      # data XDAT, NDAT
      return bin_resolution(NB,IREAD,IDUF,COMS,store,lptfile)
      """
      DER2 = vec(m_d)
      LSTART = True
      if IREAD==0:
         LSTART=False

      SMALL=1.0E-20
      COMS["FFT"].NFFT=m_d
      COMS["res_data"].N_FT = m_d//2
      xr = vec(m_d)
      yr = vec(m_d)
      er = vec(m_d)
      # copy resolution data
      xr.copy(COMS["Res"].xres.output_range(1,NB))
      yr.copy(COMS["Res"].yres.output_range(1,NB))
      er.copy(COMS["Res"].eres.output_range(1,NB))

      # rebin the resolution data as it is not on an evenly spaced grid
      XB, YB = rebin(xr,yr,er,NB,COMS["Res"].nrbin)
        
      XDMIN=COMS["DATA"].XDAT(1)
      XDMAX=COMS["DATA"].XDAT(COMS["DATA"].NDAT)
      YMAX=0.0
      YSUM=0.0
      # get indicies for valid x range
      first_index = np.argmax(XB.output() >=XDMIN)+1
      last_index = np.argmin(XB.output()<=XDMAX)+1

      # get total and max Y binned values within valid x range
      YSUM = np.sum(YB.output_range(first_index, last_index))
      YMAX = np.max(YB.output_range(first_index, last_index))

      if YSUM<SMALL:
       store.open(53,lptfile)
       store.write(53,' Ysum is too small')
       IDUF=1
       store.close(unit=53)
       return

      XBMIN, XBMAX, YB = rm_BG(XB,YB,NB,YMAX,LSTART, COMS,store, lptfile) # subtracts BG off YB
      # populate FRES with spline of evenly spaced binned data -> data to FFT later
      func=SPLINE(XB,YB,NB,0.0,0.0,DER2)
      TWOPIN=2.0*np.pi/float(COMS["FFT"].NFFT)
      COMS["FFT"].FRES.fill(0.0, COMS["FFT"].NFFT)
      XX=0.0
      DXJ=COMS["FFT"].XJ(2)-COMS["FFT"].XJ(1) # bin width
      COMS["FFT"].FRES.set(1,SPLINT(XX,func))
      SUM=COMS["FFT"].FRES(1)
      for I in get_range(1,int(COMS["FFT"].NFFT/2)):#use symmetry to have the range we need to loop
        XX=XX+DXJ
        if XX < XBMAX:
           COMS["FFT"].FRES.set(I+1,SPLINT(XX,func))
        if -XX > XBMIN:
           COMS["FFT"].FRES.set(COMS["FFT"].NFFT+1-I,SPLINT(-XX,func))
        SUM+=COMS["FFT"].FRES(I+1)+COMS["FFT"].FRES(COMS["FFT"].NFFT+1-I)
        # set phases
        COMS["FFT"].TWOPIK.set(I,TWOPIN*float(I-1)) # looks to be the phase
        COMS["res_data"].phases.set(I,TWOPIN*float(I-1))
      COMS["FFT"].TWOPIK.set(int(COMS["FFT"].NFFT/2)+1,TWOPIN*float(COMS["FFT"].NFFT/2))
      COMS["res_data"].phases.set(int(COMS["res_data"].N_FT)+1,TWOPIN*float(COMS["res_data"].N_FT))
    
      # average y value in each bin -> integral of all data is 1
      BNORM=1./(SUM*float(COMS["FFT"].NFFT))
      tmp = COMS["FFT"].FRES.output_range(1,COMS["FFT"].NFFT)
      tmp = tmp*BNORM # scale the splined data
      COMS["FFT"].FRES.copy(tmp)
      out = FOUR2(COMS["FFT"].FRES, COMS["FFT"].NFFT,1,1,0)
      COMS["FFT"].FRES.copy(flatten(out))
      COMS["res_data"].FTY.copy(out)

      # some rotations?
      for I in get_range(3,COMS["FFT"].NFFT,4):
        COMS["FFT"].FRES.set(I,-COMS["FFT"].FRES(I))
        COMS["FFT"].FRES.set(I+1, -COMS["FFT"].FRES(I+1))

      for k in get_range(2, COMS["res_data"].N_FT, 2):
          COMS["res_data"].FTY.set(k, -COMS["res_data"].FTY(k))

      # IFT of resolution
      if not LSTART:
        tmp = COMS["FFT"].FRES.output_range(end=COMS["FFT"].NFFT+2)
        COMS["FFT"].FWRK.copy(tmp)
        out = FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
        COMS["FFT"].FWRK.copy(flatten(out[0:m_d2]))


        COMS["res_data"].IFTY.copy(COMS["res_data"].FTY.output_range(end=COMS["res_data"].N_FT))
        out2 = FOUR2(COMS["res_data"].IFTY, COMS["res_data"].N_FT, 1, -1,-1)
        COMS["res_data"].IFTY.copy(out[0:m_d2])

      LSTART= True
      return XB, YB
      """

def CCHI(V,COMS, o_bgd, o_w1):
      CHI=0.0
      B1=COMS["SCL"].BSCL*V(1) # BG 1
      B2=COMS["SCL"].BSCL*V(2) # BG 2
      A0=COMS["SCL"].ASCL*V(3) # elastic peak amplitude
      DELTAX=V(4) # zero offset
      NFT2=int(COMS["FFT"].NFFT/2+1)

      # get resolution values
      fres = COMS["res_data"].FTY.output_range(end=NFT2)
      twopik = COMS["res_data"].phases.output_range(end=NFT2)
      RKEXP, RKEXP2 = CXSHFT(fres, DELTAX, twopik)#,FR2PIK,FR2PIK(1,2),NFT2)
      COMS["GRD"].FR2PIK.copy(flatten(RKEXP))
      COMS["GRD"].FR2PIK.copy(flatten(RKEXP2), 1,2)
      COMS["WORK"].WORK.fill(A0, NFT2+1)
      XNSCL=-np.log(1.0E-7)/(COMS["FFT"].TWOPIK(2)-COMS["FFT"].TWOPIK(1))
      for J in get_range(1,COMS["FIT"].NFEW):
        AJ=COMS["SCL"].ASCL*V(3+J+J)
        SIGJ=COMS["SCL"].WSCL*V(4+J+J)/COMS["SCL"].GSCL
        COMS["FIT"].EXPF.fill(0.0, NFT2, 1,J)# this might need a plus 1 to N later
        NXFT2=1+NINT(XNSCL/(np.abs(SIGJ)+1.0E-10))
        if NXFT2 > NFT2:
           NXFT2=NFT2
        for I in get_range(1,NXFT2):
          EXPIJ=np.exp(-COMS["FFT"].TWOPIK(I)*SIGJ)
          COMS["FIT"].EXPF.set(I,J, EXPIJ)          
          #if I == NXFT2:
          #print("hi", EXPIJ, COMS["FIT"].EXPF(I,J), COMS["FFT"].TWOPIK(I), SIGJ, I)

          COMS["WORK"].WORK.set(I,1, COMS["WORK"].WORK(I,1)+AJ*EXPIJ)

      tmp = VMLTRC(COMS["WORK"].WORK.output_range(1,1,end=NFT2+1),RKEXP)#,NFT2,FWRK)
      COMS["FFT"].FWRK.copy(flatten(tmp))
      tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["FFT"].FWRK.copy(flatten(tmp))
      COMS["FIT"].FIT.copy(DEGRID(COMS["FFT"].FWRK,COMS))
      X1=COMS["DATA"].XDAT(1)
      BNRM = 0.0
      if o_bgd==2:
         BNRM=(B2-B1)/(COMS["DATA"].XDAT(COMS["DATA"].NDAT)-X1)
      #avoid conflict BNORM with ModPars
      fit = COMS["FIT"].FIT.output_range(end=COMS["DATA"].NDAT)
      xdat = COMS["DATA"].XDAT.output_range(end=COMS["DATA"].NDAT)
      dat = COMS["DATA"].DAT.output_range(end=COMS["DATA"].NDAT)
      sig = COMS["DATA"].SIG.output_range(end=COMS["DATA"].NDAT)
      fit += B1+BNRM*(xdat-X1)
      diff = fit - dat
      resid = diff*sig
      CHI = np.sum(diff*resid)

      COMS["FIT"].FIT.copy(fit) # does not match at later indicies -> probably fine
      COMS["FIT"].RESID.copy(resid)

      if o_w1 == 1 and COMS["FIT"].NFEW>=1:
         RESW1D=(COMS["SCL"].WSCL*V(6)-COMS["QW1"].QW1(COMS["QW1"].ISPEC))/COMS["QW1"].SIGQW1(COMS["QW1"].ISPEC)
         CHI=CHI+2.0*pow(RESW1D,2)
      return CHI/(2.0*float(COMS["DATA"].NDAT))



def REFINE(COMS, GRAD,HESS,NP,DETLOG,INDX,COVAR,STEPSZ, o_bgd, o_w1,o_el, prog):
      # FFT FWRK, WORK WORK, GRD DDDPAR
      NFT2=int(COMS["FFT"].NFFT/2)+1
      CNORM=CCHI(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
      HESS = matrix_2(NP,NP)
      tmp = VMLTRC(COMS["WORK"].WORK.output_range(1,1,end=NFT2+1),compress(COMS["GRD"].FR2PIK.output_range(1,2,end=2*NFT2+2)))
      COMS["FFT"].FWRK.copy(flatten(tmp))
      tmp = VMLTRC(COMS["FFT"].TWOPIK.output_range(end=NFT2+1),compress(COMS["FFT"].FWRK.output_range(end=2*(NFT2+2))))
      COMS["WORK"].WORK.copy(flatten(tmp))
      
      tmp = compress(COMS["WORK"].WORK.output_range(1,1,2*(NFT2+1)))
      COMS["WORK"].WORK.copy(flatten(VMLTIC(tmp)))
      tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["FFT"].FWRK.copy(flatten(tmp))    
      COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,4)
      tmp=FOUR2(COMS["WORK"].WORK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["WORK"].WORK.copy(flatten(tmp))

      COMS["FFT"].FWRK.copy(DEGRID(COMS["WORK"].WORK.output_as_vec(),COMS))

      HESS.set(4,4,VRDOTR(COMS["FIT"].RESID.output_range(end=COMS["DATA"].NDAT+1),COMS["FFT"].FWRK.output_range(end=COMS["DATA"].NDAT+1),COMS["DATA"].NDAT-1))
      COMS["FFT"].FWRK.copy(COMS["GRD"].FR2PIK.output_range(1,1,end=COMS["FFT"].NFFT+2))
      tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["FFT"].FWRK.copy(flatten(tmp))    
      COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,3)

      COMS["FFT"].FWRK.copy(COMS["GRD"].FR2PIK.output_range(1,2,end=COMS["FFT"].NFFT+2))
      tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["FFT"].FWRK.copy(flatten(tmp))    
      COMS["WORK"].WORK.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,1)
      HESS.set(3,4,VRDOTR(COMS["FIT"].RESID.output_range(end=COMS["DATA"].NDAT+1),COMS["WORK"].WORK.output_range(1,1,end=COMS["DATA"].NDAT+2),COMS["DATA"].NDAT-1))
      HESS.set(4,3, HESS(3,4))
      for I in get_range(1,COMS["FIT"].NFEW): # something strange is happening in here...
        J=3+I+I
        AJ=COMS["FIT"].FITP(J)*COMS["SCL"].ASCL

        tmp= VMLTRC(COMS["FIT"].EXPF.output_range(1,I, end = NFT2+2),compress(COMS["GRD"].FR2PIK.output_range(1,1,2*(2+NFT2))))
        COMS["WORK"].WORK.copy(flatten(tmp))
        COMS["FFT"].FWRK.copy(COMS["WORK"].WORK.output_range(1,1,COMS["FFT"].NFFT+2))
        tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
        COMS["FFT"].FWRK.copy(flatten(tmp))   
        #COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,J)

        #tmp= VMLTRC(COMS["FFT"].TWOPIK.output_range(end = NFT2-1),compress(COMS["WORK"].WORK.output_range(1,1,2*(NFT2))))
        #COMS["FFT"].FWRK.copy(flatten(tmp))
        # up to here -> slow down is either four and/or degrid
        #COMS["WORK"].WORK.copy(COMS["FFT"].FWRK.output_range(end=COMS["FFT"].NFFT+3))
        #tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
        #COMS["FFT"].FWRK.copy(flatten(tmp))   
        #COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,J+1)# -> this messes up hessian! 
        #tmp, HESS = HESS0(HESS,COMS["FIT"].RESID.output_range(end=COMS["DATA"].NDAT),COMS["GRD"].DDDPAR.output_range(1,J+1,COMS["DATA"].NDAT+1),AJ,J) # the return vals for tmp are causing a huge slow down
        #COMS["GRD"].DDDPAR.copy(tmp, 1,J+1)
        #CALL VMLTIC(WORK,NFT2,WORK)
        #CALL VCOPY(WORK,FWRK,NFFT+2)
        #CALL FOUR2(FWRK,NFFT,1,-1,-1)
        #CALL DEGRID(FWRK,WORK(1,2))
        #CALL VRDOTR(RESID,WORK(1,2),NDAT,HESS(4,J))
        #HESS(J,4)=HESS(4,J)
        #CALL VMLTRC(TWOPIK,WORK,NFT2,FWRK)
        #CALL VCOPY(FWRK,WORK,NFFT+2)
        #CALL FOUR2(FWRK,NFFT,1,-1,-1)
        #CALL DEGRID(FWRK,WORK(1,2))
        #CALL VRDOTR(RESID,WORK(1,2),NDAT,SM)
        #HESS(4,J+1)=-AJ*SM
        #HESS(J+1,4)=HESS(4,J+1)
        #CALL VMLTIC(WORK,NFT2,WORK)
        #CALL FOUR2(WORK,NFFT,1,-1,-1)
        #CALL DEGRID(WORK,FWRK)
        #CALL VRDOTR(RESID,FWRK,NDAT,SM)
        #HESS(J+1,J+1)=-AJ*SM
      GRAD.copy(GRADPR(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS, col=2)) # Grad has different sign in index 3, I think it is "zero" (e-4)
      HESS=HESS1(NP,COMS["SCL"].SCLVEC.output_col(2),STEPSZ,COMS["FIT"].NFEW, prog, COMS,o_el, HESS=HESS) 
      if o_w1==1 and NP>6:
          #print("mksfdksfdksfdaksfd") # not tested in here
          DIF=COMS["SCL"].WSCL*COMS["FIT"].FITP(6)-COMS["QW1"].QW1(COMS["QW1"].ISPEC)
          SIG2=2.0/pow(COMS["QW"].SIGQW1(COMS["QW"].ISPEC),2)
          GRAD.set(6, GRAD(6)+SIG2*DIF*COMS["SCL"].WSCL)
          HESS.set(6,6, HESS(6,6)+SIG2*pow(COMS["SCL"].WSCL,2))
      covar_default = 1.
      if prog=='s':
          cov=2.0
      HESS, COVAR, DETLOG = INVERT(NP,INDX,covar_default, HESS)
      return HESS, COVAR, DETLOG


#***<see the fit>*******************************************************
def SEEFIT(COMS, SIGPAR,CNORM, store, lptfile):
      store.open(53, lptfile)
      ERRSCL=sqrt(CNORM)
      PRMSV = []
      SIGSV = []

      store.write(53, f' Best-fit assuming no. of quasi-elastic lines = { COMS["FIT"].NFEW}')
      store.write(53,f' >>> Normalised Chi-squared = {CNORM:.4f}')
      store.write(53, f' Background(Xmin) = {COMS["FIT"].FITP(1)*COMS["SCL"].BSCL:.3e} +- {SIGPAR(1)*COMS["SCL"].BSCL*ERRSCL:.3e}')
      store.write(53, f' Background(Xmax) = {COMS["FIT"].FITP(2)*COMS["SCL"].BSCL:.3e}  +- {SIGPAR(2)*COMS["SCL"].BSCL*ERRSCL:.3e}')
      store.write(53,f' Zero offset       = {COMS["FIT"].FITP(4)*COMS["SCL"].GSCL*1000.:.2e}  +- {SIGPAR(4)*COMS["SCL"].GSCL*ERRSCL*1000.:.2e}  ueV')
      store.write(53,' Elastic line')
      store.write(53, f'Amplitude  =   {COMS["FIT"].FITP(3)*COMS["SCL"].ASCL:.4e}  +- {SIGPAR(3)*COMS["SCL"].ASCL*ERRSCL:.4e} ')
      PRMSV.append(COMS["FIT"].FITP(3)*COMS["SCL"].ASCL)
      SIGSV.append(SIGPAR(3)*COMS["SCL"].ASCL*ERRSCL)
      for I in get_range(1,COMS["FIT"].NFEW):
        J=4+I+I
        store.write(53,f' Quasi-elastic line {I}')
        store.write(53, f' FWHM        =   {2000*COMS["FIT"].FITP(J)*COMS["SCL"].WSCL:.2e}  +- {2000*SIGPAR(J)*COMS["SCL"].WSCL*ERRSCL:.2e} ueV')
        store.write(53, f' Amplitude  =   {COMS["FIT"].FITP(J-1)*COMS["SCL"].ASCL:.4e}  +-  {SIGPAR(J-1)*COMS["SCL"].ASCL*ERRSCL:.4e}') # this seems to be wrong for just the FITP(J-1)
        PRMSV.append(COMS["FIT"].FITP(J-1)*COMS["SCL"].ASCL)
        SIGSV.append(SIGPAR(J-1)*COMS["SCL"].ASCL*ERRSCL)
        PRMSV.append(2.0*COMS["FIT"].FITP(J)*COMS["SCL"].WSCL)
        SIGSV.append(2.0*SIGPAR(J)*COMS["SCL"].WSCL*ERRSCL)
      
      store.close(53)
      return np.asarray(PRMSV), np.asarray(SIGSV)

def OUTPRM(P,C,NP,NFEW,CNORM, store, files):
      # p = param, C = covar
      #print("outfsda", NFEW) # still needs testing
      if NFEW<1 or NFEW>len(files):
          return
      for k in range(len(files)):
          store.open(k+1, files[k])
      store.write(NFEW, f'{P(3)}   {P(1)}   {P(2)}   {P(4)}')
      
      for I in get_range(5,NP,2):
        store.write(NFEW,f'{P(I)}   {P(I+1)}')
      
      CSCALE=2.0*CNORM
      for J in get_range(1,NP):
        for I in get_range(1,NP):
          C.set(I,J,CSCALE*C(I,J))
      store.write(NFEW, f'{C(3,3)}')
      store.write(NFEW, f'{C(3,5)}   {C(5,5)}')
      store.write(NFEW,f'{C(3,6)}   {C(5,6)}    {C(6,6)}')
      if NFEW > 1:
        store.write(NFEW, f'{C(3,7)}   {C(5,7)}   {C(6,7)}    {C(7,7)}')
        store.write(NFEW, f'{C(3,8)}   {C(5,8)}   {C(6,8)}   {C(7,8)}   {C(8,8)}')
      
      elif NFEW>2:
        store.write(NFEW, f'{C(3,9)}   {C(5,9)}   {C(6,9)}   {C(7,9)}   {C(8,9)}   {C(9,9)}')
        store.write(NFEW, f'{C(3,10)}   {C(5,10)}   {C(6,10)}    {C(7,10)}     {C(8,10)}   {C(9,10)}    {C(10,10)}')
      
      store.write(NFEW,' -------------------------------------------------')
      store.close(unit=1)
      store.close(unit=2)
      store.close(unit=3)

      
#**<search for one more & refine amplitudes>****************************
def SEARCH(COMS, GRAD,HESS,DPAR, INDX,COVAR, o_w1, prog, o_bgd, o_el, store, lptfile, DETLOG, Chi_func):
    #FIt FITP, HESS, COVAR, DPAR,  FFT FWRK, GRD DDDPAR
      if o_w1 == 1 and COMS["FIT"].NFEW>=1:
          COMS["FIT"].FITP.set(5, 0.1)
          COMS["FIT"].FITP.set(6, COMS["QW1"].QW1(COMS["QW1"].ISPEC)/COMS["SCL"].WSCL)
          if COMS["FIT"].NFEW == 1:
              HESS, COVAR, DPAR=REFINA(GRAD,3+COMS["FIT"].NFEW,DETLOG,INDX,COVAR, COMS, Chi_func, prog, o_bgd,o_w1, o_el, store, lptfile)
              return HESS, COVAR, DPAR
      J=4+2*COMS["FIT"].NFEW
      DXLOG=0.85
      NSRCH=NINT(np.log(5.0*COMS["SCL"].GSCL/COMS["SCL"].WSCL)/np.log(DXLOG))
      CMIN=1.0E20
      COMS["FIT"].FITP.set(J-1, 0.1)
      COMS["FIT"].FITP.set(J, 1.0)
      SIGJ=COMS["FIT"].FITP(J)
      print("moooo", NSRCH)
      for I in get_range(1,NSRCH):
        HESS, COVAR, DPAR=REFINA(GRAD,3+COMS["FIT"].NFEW,DETLOG,INDX,COVAR, COMS, Chi_func, prog, o_bgd,o_w1, o_el, store, lptfile)
        CNORM=Chi_func(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
        if CNORM<CMIN:
          CMIN=CNORM
          SIGJ=COMS["FIT"].FITP(J)
       
        COMS["FIT"].FITP.set(J, COMS["FIT"].FITP(J)*DXLOG)
#      print("SIGJ", SIGJ, CMIN, CNORM, J)
      COMS["FIT"].FITP.set(J, SIGJ)
      return REFINA(GRAD,3+COMS["FIT"].NFEW,DETLOG,INDX,COVAR, COMS, Chi_func, prog, o_bgd,o_w1, o_el, store, lptfile) 
