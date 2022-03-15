# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from quasielasticbayes.python.fortran_python import *
#from quasielasticbayes.python.constants import *
#from quasielasticbayes.python.four import *
from math import pi
import numpy as np
from quasielasticbayes.python.four import *



def DPINIT(COMS):
      I=1
      COMS["SCL"].GSCL=(COMS["FFT"].XJ(COMS["FFT"].NFFT)-COMS["FFT"].XJ(1))/float(COMS["FFT"].NFFT-1) # average XJ interval

      def data_check(XJ, XDAT, J):
          return XJ(J) - XDAT >=-5.e-6 # given a tol due to small differences
      for K in get_range(1,COMS["DATA"].NDAT):
        J = find_index((COMS["FFT"].XJ, COMS["DATA"].XDAT(K)),I,COMS["FFT"].NFFT, data_check)
        COMS["Dintrp"].IPDAT.set(K, J-1)
        # normalised XDAT-JDAT
        COMS["Dintrp"].XPDAT.set(K, (COMS["DATA"].XDAT(K)-COMS["FFT"].XJ(J-1))/COMS["SCL"].GSCL)
        I=J

def GDINIT(COMS):
      X1=COMS["DATA"].XDAT(1)
      XN=COMS["DATA"].XDAT(COMS["DATA"].NDAT)
      GNORM=1.0/(XN-X1)
      for I in get_range(1,COMS["DATA"].NDAT):
        XGNORM=(COMS["DATA"].XDAT(I)-X1)*GNORM # fraction of x range
        COMS["GRD"].DDDPAR.set(I,1, 1.0-XGNORM) # normalized offset of xdat
        COMS["GRD"].DDDPAR.set(I,2, XGNORM)

"""
***<read in the data>**************************************************
"""

def DATIN1(COMS, store, lptfile):
      SMALL=1.0E-20
      if abs(COMS["Params"].RSCL-1.0) > 0.01:
       store.open(53,lptfile)
       store.write(53,f' DATIN1; Data error-bars multiplied by: {COMS["Params"].RSCL}')
       store.close(unit=53)
      COMS["Params"].RSCL =pow(COMS["Params"].RSCL, 2)
      N=0
      # lopp over bins
      for I in get_range(COMS["Params"].IMIN,COMS["Params"].IMAX,COMS["Params"].NBIN):
       N=N+1
       XXD=0.0
       DD=0.0
       EE=0.0
       K=0
       # looping across each bin
       for J in get_range(0,COMS["Params"].NBIN-1):
        XXD=XXD+COMS["DATA"].XDAT(I+J)
        if COMS["DATA"].SIG(I+J) > SMALL:
         K=K+1
         DD=DD+COMS["DATA"].DAT(I+J)
         EE=EE+COMS["DATA"].SIG(I+J)
       # store the scaled value for the new bin
       COMS["DATA"].XDAT.set(N, COMS["Params"].BNORM*XXD)
       if K > 0:
        # if large sigma(s) exist in bin assume it dominates -> throw away everyting else
        COMS["DATA"].DAT.set(N, COMS["Params"].BNORM*DD)
        COMS["DATA"].SIG.set(N, 2.0*float(K*K)/(EE*COMS["Params"].RSCL))
       else:
        COMS["DATA"].DAT.set(N, 0.0)
        COMS["DATA"].SIG.set(N, 0.0)

      COMS["DATA"].NDAT=N

def DATIN(IREAD,DTNORM, efix, ntc, COMS, store,lptfile):
      IDUF=0
      SMALL=1.0E-10
      DSUM=0.0
     
      COS2TH=2.0*np.cos(COMS["DATA"].theta(IREAD)*pi/180.0)
      QQ=efix+efix-COS2TH*abs(efix)
      COMS["DATA"].QAVRG.set(IREAD, 0.69469*np.sqrt(QQ))
      store.open(53,lptfile)
      store.write(53,' ----------------------------------------------------')
      store.write(53, f' Group {IREAD},  theta = {COMS["DATA"].theta(IREAD)},  Q = {COMS["DATA"].QAVRG(IREAD)}')
      store.write(53,' ----------------------------------------------------')
      DTNRM=DTNORM(IREAD)
      NDAT=ntc-1

      COMS["DATA"].XDAT.copy(COMS["DATA"].xin.output_range(end=COMS["DATA"].NDAT))
      COMS["DATA"].DAT.copy(COMS["DATA"].yin.output_range(end=COMS["DATA"].NDAT))

      for I in get_range(1,COMS["DATA"].NDAT):
       COMS["DATA"].SIG.set(I, COMS["DATA"].ein(I)*DTNRM)

       if COMS["DATA"].SIG(I) > SMALL:
        COMS["DATA"].SIG.set(I, pow(COMS["DATA"].SIG(I),2))
        DSUM=DSUM+COMS["DATA"].DAT(I)
       else:
        COMS["DATA"].SIG.set(I, 0.0)

      if DSUM < SMALL:
         IDUF=1
      DATIN1(COMS, store, lptfile)
      store.close(unit=53)
      return IDUF

# dump data to a file
def FileInit(Nu, ISP, COMS, store, files):

      for n in range(Nu):
       store.open(n, files[n])
       store.write(n,f'{COMS["DATA"].QAVRG(ISP)}, {COMS["SCL"].ASCL}, {COMS["SCL"].WSCL}, {COMS["SCL"].BSCL}, {COMS["SCL"].GSCL}')
       store.close(unit=n)

def PRINIT(NQMAX,IXSCAL,COMS,store, prog,lptfile, o_bgd ):
      """
      This seems to assume that a large error -> its BG measurment
      If its a peak the data should dominate the measurment and the error is small
      """

      if IXSCAL <= 1:
         SMALL=1.0E-10
         SM=0.0
         NSUM=0
         # get upto 10 "large" values in SIG from first 20 values
         for I in get_range(1,20):
          if COMS["DATA"].SIG(I) >=SMALL:
           NSUM=NSUM+1
           SM=SM+abs(COMS["DATA"].DAT(I))
          if NSUM >= 10:
              break
         # assume large values dominate and normalise
         BSCL1=SM/float(NSUM)

         SM=0.0
         NSUM=0
         # get upto 10 "large" values in SIG from last 20 values
         for I in get_range(1,20):
          if COMS["DATA"].SIG(COMS["DATA"].NDAT-I+1) >= SMALL:
           NSUM=NSUM+1
           SM=SM+abs(COMS["DATA"].DAT(COMS["DATA"].NDAT-I+1))
          
          if NSUM >= 10:
             break

         COMS["SCL"].BSCL=SM/float(NSUM)
         # use the smaller normalization
         if BSCL1<COMS["SCL"].BSCL:
            COMS["SCL"].BSCL=BSCL1
         COMS["SCL"].BSCL=COMS["SCL"].BSCL/2.0

         if o_bgd == 0:
             COMS["SCL"].BSCL=0.0 # zero background

         # get the mod largest start/end value
         AXMAX=abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))
         if abs(COMS["DATA"].XDAT(1))>AXMAX:
            AXMAX=abs(COMS["DATA"].XDAT(1))

         COMS["SCL"].WSCL=AXMAX/3.0 # no idea why 3

         MK=0
         SM=0.0
         SUMSIG=0.0
         for I in get_range(1,COMS["DATA"].NDAT-1):
          # if not BG
          if COMS["DATA"].SIG(I) >= SMALL:
           MK=MK+1 
           # remove BG and scale by bin width
           SM=SM+(COMS["DATA"].DAT(I)-COMS["SCL"].BSCL)*(COMS["DATA"].XDAT(I+1)-COMS["DATA"].XDAT(I))
           SUMSIG=SUMSIG+np.sqrt((2.0/COMS["DATA"].SIG(I))) # s = sigma sqrt{\frac{2}{sigma(I)} }

         # scale the average of the non-BG data by total bins
         COMS["SCL"].ASCL=SM*float(COMS["DATA"].NDAT)/float(MK)
         # average error of non-BG data
         SUMSIG=SUMSIG/float(MK)
         # scale av error (not sure why sqrt)
         SUMSIG=(COMS["DATA"].XDAT(COMS["DATA"].NDAT)-COMS["DATA"].XDAT(1))*SUMSIG/np.sqrt(float(MK))
         store.open(53,lptfile)
         # if measurments are less then errors
         if COMS["SCL"].ASCL < SUMSIG:
           store.write(53,' qlm> *** Estimate of Amax is being set to lower bound!')
           store.write(53,f' ( {COMS["SCL"].ASCL} --> {SUMSIG} )')
           COMS["SCL"].ASCL=SUMSIG
         store.write(53,' ----------------------------------------------------')
         store.close(unit=53)
         COMS["SCL"].ASCL=COMS["SCL"].ASCL/COMS["SCL"].GSCL # rescale 
         # no idea - fit parameters?
         COMS["FIT"].FITP.set(1, 1.0)
         COMS["FIT"].FITP.set(2, 1.0)
         COMS["FIT"].FITP.set(3, 0.5)
         # seems to be storing the max values and BG for latter
         for I in get_range(1,2):
           COMS["SCL"].SCLVEC.set(I,1,COMS["SCL"].BSCL)
           COMS["SCL"].SCLVEC.set(I,2,COMS["SCL"].BSCL)
         COMS["SCL"].SCLVEC.set(3,1,COMS["SCL"].ASCL)
         COMS["SCL"].SCLVEC.set(3,2,COMS["SCL"].ASCL)
      COMS["SCL"].SCLVEC.set(4,2,1.0)
      if prog == 'w':
       COMS["SCL"].SCLVEC.set(4,1, COMS["SCL"].ASCL)
       COMS["SCL"].SCLVEC.set(5,2, COMS["SCL"].ASCL)
       COMS["SCL"].SCLVEC.set(6,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
       COMS["SCL"].SCLVEC.set(7,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
      else:
       print("hi", NQMAX)
       for I in get_range(1,NQMAX):
        COMS["SCL"].SCLVEC.set(3+I,1, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(3+I+I,2, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(4+I+I,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
      COMS["FIT"].NFEW=0
