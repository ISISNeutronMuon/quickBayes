# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from quasielasticbayes.python.fortran_python import *
#from quasielasticbayes.python.constants import *
#from quasielasticbayes.python.four import *
from math import pi, log10, sqrt
import numpy as np
from quasielasticbayes.python.four import *
from quasielasticbayes.python.util import *



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


def CXSHFT(RK,DX,TWOPIK):
    XX = TWOPIK*DX
    XC = np.cos(XX)+ 1j*np.sin(XX)
    RKEXP = RK*XC
    RKEXP2 = VMLTRC(TWOPIK,RKEXP)
    RKEXP2=VMLTIC(RKEXP2)
    return RKEXP, RKEXP2

def VMLTRC(R,C):
    A = R*C.real
    B = R*C.imag
    C = A + 1j*B
    return A + 1j*B
#C     -------------------------
def VMLTIC(C):
    return C*1j

# seems to rescale the YGRD data
def DEGRID(YGRD, COMS):
      YDAT = []
      for I in get_range(1,COMS["DATA"].NDAT):
        J=int(COMS["Dintrp"].IPDAT(I))
        YDAT.append(YGRD(J)+COMS["Dintrp"].XPDAT(I)*(YGRD(J+1)-YGRD(J)))
      return np.asarray(YDAT)

def VRDOTR(A,B):
      return np.sum(A*B)

# weights the SCLVEC
def GRADPR(RESID,NDAT,NP,SCLVEC, COMS,col=1):
      GRAD = []
      for I in get_range(1,NP):
        SM=VRDOTR(RESID.output_range(end=NDAT),COMS["GRD"].DDDPAR.output_range(1,I, end=NDAT+1))#,NDAT,SM)
        GRAD.append(SCLVEC(I,col)*SM)
      return np.asarray(GRAD)


def HESS0(HESS, RESID, DDDPAR,AJ,J):
      SM = np.sum(-RESID*DDDPAR)
      HESS.set(J,J+1, SM)
      HESS(J+1,J,SM) # symmetric matrix
      return -DDDPAR*AJ

# if HESS is None, then create an NP by NP Hessian matrix
def HESS1(NP,SCLvec,STEPSZ,NFEW, prog, COMS, HESS=None):
      if HESS is None:
          HESS = matrix_2(NP,NP)
      for J in get_range(1,NP):
        for I in get_range(J,NP):
          SM=np.sum(COMS["DATA"].SIG.output_range(end=COMS["DATA"].NDAT)*COMS["GRD"].DDDPAR.output_range(1,I,COMS["DATA"].NDAT+1)*COMS["GRD"].DDDPAR.output_range(1,J,COMS["DATA"].NDAT+1))
          HESS.set(I,J, (HESS(I,J)+SM)*SCLvec[I-1]*SCLvec[J-1])
          HESS.set(J,I, HESS(I,J)) # symmetric hessian

      BEEFUP=2.0/(STEPSZ*STEPSZ)
      for I in get_range(1,NP):
        HESS.set(I,I, HESS(I,I)+BEEFUP)

      if prog=='l' or prog=='s':
       if NFEW>0: # option for elastic peak
        if o_el==0:
           HESS.set(3,3,2.0E8)
      return HESS

def INVERT(NP, INDX, covar_default, HESS=None, COVAR=None):
    if HESS is None:
        HESS = matrix_2(NP,NP)
    if COVAR is None:
        COVAR = matrix_2(NP,NP)
    SMALL=1.E-20
    DETLOG=0.0
    COVAR.fill(0.0, NP*NP)
    for I in get_range(1,NP):
        COVAR.set(I,I, covar_default)
    INDX,D=LUDCMP(HESS,NP,NP)
    for I in get_range(1,NP):
      DETLOG=DETLOG+log10(abs(HESS(I,I))+SMALL)
        
    for I in get_range(1,NP):
        tmp = LUBKSB(HESS,NP,NP,INDX,COVAR.output_from(1,I))
        COVAR.copy(tmp, 1,I)
    return HESS, COVAR, DETLOG

def MLTMXV(P,OP,N):
      D = vec(N*N) # could be N long?
      for K in get_range(1,N):
        SM=0.#np.sum(p_vec)
        for J in get_range(1,N):
          SM=SM+OP(J,K)*P(J)
        #end do
        D.set(K,SM)
      return D


def NEWEST(COVAR,GRAD,NP,NFEW,FITP, prog, store,lptfile):

      if prog=='w':
       mp=3
      else:
       mp=2
      
      DPAR = MLTMXV(GRAD,COVAR,NP)# determinenet of grad and covar
      if NP == 4+mp*NFEW:
        for I in get_range(1,NP):
          FITP.set(I, FITP(I)-DPAR(I))
    
      elif NP == 3+NFEW:
        for I in get_range(1,3):
          FITP.set(I,FITP(I)-DPAR(I))
      
        for I in get_range(1,NFEW):
          J=I+3
          FITP.set(J+I, FITP(J+I)-DPAR(J))
      
      else:
         store.open(53,lptfile )
         store.write(53,' NEWEST : Something wrong here folks!')
         store.close(53)
      return DPAR

# should we define grad in here too?
def REFINA(GRAD,NP,DETLOG,INDX,COVAR, COMS, CNORM_FUNC, prog, o_bgd,o_w1, o_el, store, lptfile):
      NFT2=COMS["FFT"].NFFT/2+1
      CNORM=CNORM_FUNC(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
      #HESS.fill(0.0,NP*NP)
      COMS["FFT"].FWRK.copy(COMS["GRD"].FR2PIK.output_range(1,1,COMS["FFT"].NFFT+2))
      tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      COMS["FFT"].FWRK.copy(flatten(tmp))
      COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,3)
      for I in get_range(1,COMS["FIT"].NFEW):
        tmp = VMLTRC(COMS["FIT"].EXPF.output_range(1,I,end=NFT2+1),COMS["GRD"].FR2PIK.output_range(1,1,end=NFT2+1))#,NFT2,FWRK)
        COMS["FFT"].FWRK.copy(flatten(tmp))
        tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
        COMS["FFT"].FWRK.copy(flatten(tmp))
        COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,3+I)
      GRAD.copy(GRADPR(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS))
      HESS=HESS1(NP,COMS["SCL"].SCLVEC.output(),0.3,COMS["FIT"].NFEW, prog, COMS, HESS=None) # create HESS matrix
      covar_default = 1
      if prog == 's':
          covar_default = 2
      HESS, COVAR, DETLOG = INVERT(NP,INDX,covar_default, HESS)
      # changes FITP
      DPAR = NEWEST(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)
      CNORM=CNORM_FUNC(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
      GRAD.copy(GRADPR(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS))
      DPAR = NEWEST(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)

      return HESS, COVAR, DPAR

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
       for I in get_range(1,NQMAX):
        COMS["SCL"].SCLVEC.set(3+I,1, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(3+I+I,2, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(4+I+I,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
      COMS["FIT"].NFEW=0

def ERRBAR(COVAR,NP):
    # gets the error bars from the diag of the covarience matrix
    SMALL=1.0E-20
    SIGPAR = vec(NP)
    for I in get_range(1,NP):
        SIGPAR.set(I, sqrt(2.0*abs(COVAR(I,I))+SMALL))
    return SIGPAR
