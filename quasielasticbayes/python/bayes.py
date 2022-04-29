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

def bin_offsets(COMS):
      I=1
      # If NFFT is large the minus one makes little difference and this is roughly the average bin width -> but bin width is conceptually easier
      COMS["SCL"].GSCL=(COMS["FFT"].XJ(COMS["FFT"].NFFT)-COMS["FFT"].XJ(1))/float(COMS["FFT"].NFFT-1)#  average bin width

      def data_check(XJ, XDAT, J):
          return XJ(J) - XDAT >=-5.e-6 # given a tol due to small differences
      for K in get_range(1,COMS["DATA"].NDAT):
        # get first XJ that lines up with bin at XDAT(K)
        J = find_index((COMS["FFT"].XJ, COMS["DATA"].XDAT(K)),I,COMS["FFT"].NFFT, data_check)
        COMS["Dintrp"].IPDAT.set(K, J-1) # record how new bins line up with original bins
        # normalised XDAT-JDAT
        COMS["Dintrp"].XPDAT.set(K, (COMS["DATA"].XDAT(K)-COMS["FFT"].XJ(J-1))/COMS["SCL"].GSCL) # fractional offset in the bin values
        I=J


@deprecated
def DPINIT(COMS):
      #bin_offsets(COMS)
      I=1
      COMS["SCL"].GSCL=(COMS["FFT"].XJ(COMS["FFT"].NFFT)-COMS["FFT"].XJ(1))/float(COMS["FFT"].NFFT-1) # number of bins
      def data_check(XJ, XDAT, J):
          return XJ(J) - XDAT >=-5.e-6 # given a tol due to small differences
      for K in get_range(1,COMS["DATA"].NDAT):
        # get first XJ that lines up with bin at XDAT(K)
        J = find_index((COMS["FFT"].XJ, COMS["DATA"].XDAT(K)),I,COMS["FFT"].NFFT, data_check)
        COMS["Dintrp"].IPDAT.set(K, J-1) # record how new bins line up with original bins
        # normalised XDAT-JDAT
        COMS["Dintrp"].XPDAT.set(K, (COMS["DATA"].XDAT(K)-COMS["FFT"].XJ(J-1))/COMS["SCL"].GSCL) # offset in the bin values
        I=J
       

def normalize_x_range(COMS):
      X1=COMS["DATA"].XDAT(1)
      X_range=COMS["DATA"].XDAT(COMS["DATA"].NDAT) -X1
      norm_for_x_range = 1.0/(X_range)
      for I in get_range(1,COMS["DATA"].NDAT):
        normalized_x = (COMS["DATA"].XDAT(I)-X1)*norm_for_x_range # fraction of x range
        COMS["GRD"].DDDPAR.set(I,1, 1.0-normalized_x)
        COMS["GRD"].DDDPAR.set(I,2, normalized_x)


@deprecated
def GDINIT(COMS):
      #normalize_x_range(COMS)
      X1=COMS["DATA"].XDAT(1)
      XN=COMS["DATA"].XDAT(COMS["DATA"].NDAT)
      GNORM=1.0/(XN-X1)
      for I in get_range(1,COMS["DATA"].NDAT):
        XGNORM=(COMS["DATA"].XDAT(I)-X1)*GNORM # fraction of x range
        COMS["GRD"].DDDPAR.set(I,1, 1.0-XGNORM)
        COMS["GRD"].DDDPAR.set(I,2, XGNORM)
      

#fres, DELTAX, 
def complex_shift(RK,DX,TWOPIK):
    XX = TWOPIK*DX # oscillation term (bin_width*phase factor)
    XC = np.cos(XX)+ 1j*np.sin(XX)
    RKEXP = RK*XC # FT(resolution)*oscillations -> Fourier cosine plus Fourier Sin transformation args
    # this multiplies the above by 2i*pi*k -> normalisation from FT?
    RKEXP2 = TWOPIK*RKEXP*1j
    return RKEXP, RKEXP2

@deprecated
def CXSHFT(RK,DX,TWOPIK):
    XX = TWOPIK*DX # oscillation term (offset*phase factor)
    XC = np.cos(XX)+ 1j*np.sin(XX)
    RKEXP = RK*XC # FT(resolution)*oscillations
    RKEXP2 = VMLTRC(TWOPIK,RKEXP) # multiply together to get even and odd part
    RKEXP2=VMLTIC(RKEXP2) # times by i -> phase factors are missing an i
    return RKEXP, RKEXP2


def VMLTRC(R,C):
    return R*C

def VMLTIC(C):
    return C*1j

# shift bin values onto new grid
def bin_shift( y_grid, COMS): # is this the slow down?
      y_shifted = []
      for I in get_range(1,COMS["DATA"].NDAT):
        J=int(COMS["Dintrp"].IPDAT(I)) # get the index that says where the shift value is
        fractional_x_shift = COMS["Dintrp"].XPDAT(I) 
        # get fractions of original bins in the new shifted bin add sum (e.g fractional_x_shift = 0.2)
        y_shifted.append(y_grid[J-1]*(1-fractional_x_shift) + fractional_x_shift*y_grid[J])
      return np.asarray(y_shifted)


# seems to rescale the YGRD data
@deprecated
def DEGRID(YGRD, COMS): # is this the slow down?
      return bin_shift(YGRD.output(), COMS)
      #YDAT = []
      #for I in get_range(1,COMS["DATA"].NDAT):
      #  J=int(COMS["Dintrp"].IPDAT(I))
      #  YDAT.append(YGRD(J)+COMS["Dintrp"].XPDAT(I)*(YGRD(J+1)-YGRD(J)))
      #return np.asarray(YDAT)

def VRDOTR(A,B,N):
    sm = 0.0
    for j in get_range(1,N):
        sm += A[j-1]*B[j-1]
    #return np.sum(A*B)
    return sm


def construct_gradients(RESID,NDAT,NP,SCLVEC, COMS,col=1):
      GRAD = []
      # construct the gradient by product of residulas, parameter and function evaluation
      for I in get_range(1,NP):
        SM = np.sum(RESID.output_range(end=NDAT)*COMS["GRD"].DDDPAR.output_range(1,I, end=NDAT+1))#,NDAT,SM)
        GRAD.append(SCLVEC(I,col)*SM)
      return np.asarray(GRAD)

# weights the SCLVEC
#@deprecated
def GRADPR(RESID,NDAT,NP,SCLVEC, COMS,col=1):
      return construct_gradients(RESID,NDAT,NP,SCLVEC, COMS,col)
      #GRAD = []
      # construct the gradient by product of residulas, parameter and function evaluation
      #for I in get_range(1,NP):
      #  SM=VRDOTR(RESID.output_range(end=NDAT+1),COMS["GRD"].DDDPAR.output_range(1,I, end=NDAT+2), NDAT-1)#,NDAT,SM)
      #  GRAD.append(SCLVEC(I,col)*SM)
      #  #print('grad check', SM, I, NP, NDAT)
      #print()
      #return np.asarray(GRAD)


def HESS0(HESS, RESID, DDDPAR,AJ,J):
      SM = np.sum(-RESID*DDDPAR)
      HESS.set(J,J+1, SM)
      HESS.set(J+1,J,SM) # symmetric matrix
      return -DDDPAR*AJ, HESS

# if HESS is None, then create an NP by NP Hessian matrix
def make_hessian(N_params,scale_vector,step,N_QE, prog, COMS, o_el, HESS=None):
      if HESS is None:
          HESS = matrix_2(N_params ,N_params)
      for J in get_range(1,N_params):
        for I in get_range(J,N_params):
          # sum weights*evaluation for individual function_I * evaluation for individual function_J 
          SM=np.sum(COMS["DATA"].SIG.output_range(end=COMS["DATA"].NDAT)*COMS["GRD"].DDDPAR.output_range(1,I,COMS["DATA"].NDAT+1)*COMS["GRD"].DDDPAR.output_range(1,J,COMS["DATA"].NDAT+1))
          HESS.set(I,J, +SM*scale_vector[I-1]*scale_vector[J-1]) # this uses the previous value to give some influence to the history of the hessian
          HESS.set(J,I, HESS(I,J)) # symmetric hessian

      BEEFUP=2.0/(step*step)
      for I in get_range(1,N_params):
        HESS.set(I,I, HESS(I,I)+BEEFUP)

      if prog=='l' or prog=='s':
       if N_QE>0: # option for elastic peak
        if o_el==0:
           HESS.set(3,3,2.0E8)
      return HESS

#@deprecated
def HESS1(NP,SCLvec,STEPSZ,NFEW, prog, COMS,o_el, HESS=None):
      return make_hessian(NP,SCLvec,STEPSZ,NFEW, prog, COMS,o_el, HESS)
      #if HESS is None:
      #    HESS = matrix_2(NP,NP)
      #for J in get_range(1,NP):
      #  for I in get_range(J,NP):
      #    SM=np.sum(COMS["DATA"].SIG.output_range(end=COMS["DATA"].NDAT)*COMS["GRD"].DDDPAR.output_range(1,I,COMS["DATA"].NDAT+1)*COMS["GRD"].DDDPAR.output_range(1,J,COMS["DATA"].NDAT+1))
      #    HESS.set(I,J, (HESS(I,J)+SM)*SCLvec[I-1]*SCLvec[J-1])
      #    HESS.set(J,I, HESS(I,J)) # symmetric hessian
      #
      #BEEFUP=2.0/(STEPSZ*STEPSZ)
      #for I in get_range(1,NP):
      #  HESS.set(I,I, HESS(I,I)+BEEFUP)

      #if prog=='l' or prog=='s':
      # if NFEW>0: # option for elastic peak
      #  if o_el==0:
      #     HESS.set(3,3,2.0E8)
      #return HESS

def INVERT(NP, INDX, covar_default, HESS=None, COVAR=None):
    if HESS is None:
        HESS = matrix_2(NP,NP)
    if COVAR is None:
        COVAR = matrix_2(NP,NP)
    SMALL=1.E-20
    DETLOG=0.0
    COVAR.fill(0.0, NP*NP)
    # set diagonal for covariance matrix
    for I in get_range(1,NP):
        COVAR.set(I,I, covar_default)
    INDX,D, HESS=LUDCMP(HESS,NP,NP)
    for I in get_range(1,NP):
      DETLOG=DETLOG+log10(abs(HESS(I,I))+SMALL)
        
    #print("nooo", NP, INDX.output())
    for I in get_range(1,NP):
        tmp = LUBKSB(HESS,NP,NP,INDX,COVAR.output_col(I))
        COVAR.copy(tmp, 1,I)
        #for j in get_range(1,NP):
        #    print("tmp", j, tmp[j-1], len(tmp))
        
    return HESS, COVAR, DETLOG

def matrix_times_vector(grad,covar,N_params):
      result = vec(N_params*N_params)
      for K in get_range(1,N_params):
        element=np.sum(covar.output_range(1,K,end=N_params)*grad.output_range(end=N_params-1))
        result.set(K,element)
      return result

@deprecated
def MLTMXV(P,OP,N):
      return matrix_times_vector(P, OP, N)
      #D = vec(N*N) # could be N long?
      #for K in get_range(1,N):
      #  SM=0.#np.sum(p_vec)
      #  for J in get_range(1,N):
      #    SM=SM+OP(J,K)*P(J)
      #  #end do
      #  D.set(K,SM)
      #return D


def update_fit_params(COVAR,GRAD,N_params,N_EP,fit_params, prog, store,lptfile):
      # FIT FITP
      mp = 2 # number of params per elastic peak
      if prog=='w':
       mp=3
      
      d_params = matrix_times_vector(GRAD,COVAR,N_params)# determinenet of grad and covar
      # total number of parameters check
      if N_params == 4+mp*N_EP:
        for I in get_range(1,N_params):
          fit_params.set(I, fit_params(I)-d_params(I)) # adjust the fit parameters based on the grads and covar
    
      elif N_params == 3+N_EP:
        for I in get_range(1,3): # these are the BG and elastic parameters
          fit_params.set(I,fit_params(I)-d_params(I))
      
        for I in get_range(1,N_EP): # elastic peak
          J=I+3
          fit_params.set(J+I, fit_params(J+I)-d_params(J))
      
      else:
         store.open(53,lptfile )
         store.write(53,' update_fit_params (NEWEST): Something wrong here folks!')
         store.close(53)
      return d_params, fit_params.output()

@deprecated
def NEWEST(COVAR,GRAD,NP,NFEW,FITP, prog, store,lptfile):
      return update_fit_params(COVAR,GRAD,NP,NFEW,FITP, prog, store,lptfile)
      # FIT FITP
      #if prog=='w':
      # mp=3
      #else:
      # mp=2
      # 
      #DPAR = MLTMXV(GRAD,COVAR,NP)# determinenet of grad and covar
      #if NP == 4+mp*NFEW:
      #  for I in get_range(1,NP):
      #    FITP.set(I, FITP(I)-DPAR(I))
    
      #elif NP == 3+NFEW:
      #  for I in get_range(1,3):
      #    FITP.set(I,FITP(I)-DPAR(I))
      
      #  for I in get_range(1,NFEW):
      #    J=I+3
      #    FITP.set(J+I, FITP(J+I)-DPAR(J))
      
      #else:
      #   store.open(53,lptfile )
      #   store.write(53,' NEWEST : Something wrong here folks!')
      #   store.close(53)
      #return DPAR


def refine_param_values(GRAD,NP,DETLOG,INDX,COVAR, COMS, make_fit_and_chi_func, prog, o_bgd,o_w1, o_el, store, lptfile):
      # HESS, COVAR, DPAR, FFT FWRK, GRD DDDPAR, FIT FITP
      NFT2 = 1 + COMS["FFT"].NFFT//2
      CNORM = make_fit_and_chi_func(COMS["FIT"].FITP,COMS, o_bgd, o_w1)

      resolution_FT = COMS["GRD"].FR2PIK.output_range(1,1,COMS["FFT"].NFFT+2) 
      resolution_FT = np.pad(resolution_FT, pad_width=(0, m_d2 - len(resolution_FT)%m_d2),mode= 'constant')
      resolution=flatten(FOUR2_IFT(resolution_FT,COMS["FFT"].NFFT,1,-1)) # resolution*osc in origianl domain

      COMS["GRD"].DDDPAR.copy(bin_shift(resolution,COMS),1,3) # store resolution in the sample original bins

      # for each inelastic peak, convolve them with the resolution and transform back to original sampling and domain
      for I in get_range(1,COMS["FIT"].NFEW):
        convolution = COMS["FIT"].EXPF.output_range(1,I,end=NFT2+1)*compress(COMS["GRD"].FR2PIK.output_range(1,1,end=2*(NFT2+1))) #resolution * inelastic peak
        peak_in_original_domain = flatten(FOUR2_IFT(convolution,COMS["FFT"].NFFT,1,-1))
        peak_in_original_domain = bin_shift(peak_in_original_domain, COMS)
        COMS["GRD"].DDDPAR.copy(peak_in_original_domain,1,3+I) # store the individual peaks convolved with resolution func

      GRAD.copy(construct_gradients(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS))
      ################################################################################
      # up to here
      ###################################################################################   
      HESS=make_hessian(NP,COMS["SCL"].SCLVEC.output(),0.3,COMS["FIT"].NFEW, prog, COMS,o_el, HESS=None) # create HESS matrix function
      covar_default = 1
      if prog == 's':
          covar_default = 2
      HESS, COVAR, DETLOG = INVERT(NP,INDX,covar_default, HESS) # solves the covar and hessian matrix

      
      DPAR, new_params = update_fit_params(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)
      COMS["FIT"].FITP.copy(new_params)
      CNORM=make_fit_and_chi_func(COMS["FIT"].FITP,COMS, o_bgd, o_w1)

      GRAD.copy(construct_gradients(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS)) # function
      DPAR, new_params = update_fit_params(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile) # function
      COMS["FIT"].FITP.copy(new_params)

      return HESS, COVAR, DPAR, DETLOG


# should we define grad in here too?
@deprecated 
def REFINA(GRAD,NP,DETLOG,INDX,COVAR, COMS, CNORM_FUNC, prog, o_bgd,o_w1, o_el, store, lptfile):
      return refine_param_values(GRAD,NP,DETLOG,INDX,COVAR, COMS, CNORM_FUNC, prog, o_bgd,o_w1, o_el, store, lptfile)
      # HESS, COVAR, DPAR, FFT FWRK, GRD DDDPAR, FIT FITP
      #NFT2=int(COMS["FFT"].NFFT/2+1)
      #CNORM=CNORM_FUNC(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
      #HESS.fill(0.0,NP*NP)
      #COMS["FFT"].FWRK.copy(COMS["GRD"].FR2PIK.output_range(1,1,COMS["FFT"].NFFT+2))
      #tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      #COMS["FFT"].FWRK.copy(flatten(tmp))
      #COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,3)
      #for I in get_range(1,COMS["FIT"].NFEW):
      #  tmp = VMLTRC(COMS["FIT"].EXPF.output_range(1,I,end=NFT2+1),compress(COMS["GRD"].FR2PIK.output_range(1,1,end=2*(NFT2+1))))#,NFT2,FWRK)
      #  COMS["FFT"].FWRK.copy(flatten(tmp))
      #  tmp=FOUR2(COMS["FFT"].FWRK,COMS["FFT"].NFFT,1,-1,-1)
      #  COMS["FFT"].FWRK.copy(flatten(tmp))
      #  tmp = DEGRID(COMS["FFT"].FWRK,COMS)
      #  COMS["GRD"].DDDPAR.copy(DEGRID(COMS["FFT"].FWRK,COMS),1,3+I)

      #GRAD.copy(GRADPR(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS))
      #HESS=HESS1(NP,COMS["SCL"].SCLVEC.output(),0.3,COMS["FIT"].NFEW, prog, COMS,o_el, HESS=None) # create HESS matrix
      #covar_default = 1
      #if prog == 's':
      #    covar_default = 2
      #HESS, COVAR, DETLOG = INVERT(NP,INDX,covar_default, HESS)
      # changes FITP
      #DPAR = NEWEST(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)
      #CNORM=CNORM_FUNC(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
      #GRAD.copy(GRADPR(COMS["FIT"].RESID,COMS["DATA"].NDAT,NP,COMS["SCL"].SCLVEC, COMS))
      #DPAR = NEWEST(COVAR,GRAD,NP,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)

      #return HESS, COVAR, DPAR

"""
***<read in the data>**************************************************
"""


def bin_and_filter_sample_data(COMS, store, lptfile): # this causes a small difference
      SMALL=1.0E-20
      if abs(COMS["Params"].RSCL-1.0) > 0.01:
       store.open(53,lptfile)
       store.write(53,f' DATIN1; Data error-bars multiplied by: {COMS["Params"].RSCL}')
       store.close(unit=53)
      COMS["Params"].RSCL =pow(COMS["Params"].RSCL, 2)
      N=0
      # lopp over original bins
      for I in get_range(COMS["Params"].IMIN,COMS["Params"].IMAX,COMS["Params"].NBIN):
       N+=1
       bin_width = 0.0
       y_data_in_bin = 0.0
       sigma_data_in_bin = 0.0
       K=0
       # get values from rebinned data
       for J in get_range(0,COMS["Params"].NBIN-1):
        bin_width += COMS["DATA"].XDAT(I+J)
        if COMS["DATA"].SIG(I+J) > SMALL:
         K=K+1
         y_data_in_bin += COMS["DATA"].DAT(I+J)
         sigma_data_in_bin += COMS["DATA"].SIG(I+J)

       if K > 0:
        # if large sigma(s) exist in bin assume it dominates -> throw away everyting else
        y_data_in_bin = COMS["Params"].BNORM*y_data_in_bin
        sigma_data_in_bin = 2.0*float(K*K)/(sigma_data_in_bin*COMS["Params"].RSCL)
       else:
            y_data_in_bin = 0.0
            sigma_data_in_bin = 0.0
       # store the scaled value for the new bin
       COMS["DATA"].DAT.set(N, y_data_in_bin)
       COMS["DATA"].SIG.set(N, sigma_data_in_bin)
       COMS["DATA"].XDAT.set(N, COMS["Params"].BNORM*bin_width)

      COMS["DATA"].NDAT=N

@deprecated
def DATIN1(COMS, store, lptfile):
      #bin_and_filter_sample_data(COMS, store, lptfile)
      
      SMALL=1.0E-20
      if abs(COMS["Params"].RSCL-1.0) > 0.01:
       store.open(53,lptfile)
       store.write(53,f' DATIN1; Data error-bars multiplied by: {COMS["Params"].RSCL}')
       store.close(unit=53)
      COMS["Params"].RSCL =pow(COMS["Params"].RSCL, 2)
      N=0
      # lopp over original bins
      for I in get_range(COMS["Params"].IMIN,COMS["Params"].IMAX,COMS["Params"].NBIN):
       N=N+1
       XXD=0.0
       DD=0.0
       EE=0.0
       K=0
       # get values from rebinned data
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
      
      

def calculate_sample_bins(IREAD,DTNORM, efix, ntc, COMS, store,lptfile):
      IDUF=0
      SMALL=1.0E-10
      DSUM=0.0

      # get Average Q value
      COS2TH=2.0*np.cos(COMS["DATA"].theta(IREAD)*pi/180.0)
      QQ=efix+efix-COS2TH*abs(efix)
      COMS["DATA"].QAVRG.set(IREAD, 0.69469*np.sqrt(QQ))
      store.open(53,lptfile)
      store.write(53,' ----------------------------------------------------')
      store.write(53, f' Group {IREAD},  theta = {COMS["DATA"].theta(IREAD):10.5f},  Q = {COMS["DATA"].QAVRG(IREAD):10.5f}')
      store.write(53,' ----------------------------------------------------')
      store.close(unit=53)
      
      DTNRM=DTNORM(IREAD)
      COMS["DATA"].NDAT=ntc-1
      COMS["sample_data"].N=ntc-1

      # load sample data
      # these are already in sample data COM
      COMS["DATA"].XDAT.copy(COMS["DATA"].xin.output_range(end=COMS["sample_data"].N))
      COMS["DATA"].DAT.copy(COMS["DATA"].yin.output_range(end=COMS["sample_data"].N))

      for I in get_range(1,COMS["DATA"].NDAT):
       sigma = COMS["sample_data"].e_data(I)*DTNRM
       # only record values if errors are above a tol
       if sigma > SMALL:
        sigma = pow(sigma,2)
        DSUM += COMS["sample_data"].y_data(I)
       else:
        sigma = 0.0
       COMS["DATA"].SIG.set(I, sigma)
        
      # if no peak data
      if DSUM < SMALL:
         IDUF=1
      bin_and_filter_sample_data(COMS, store, lptfile)

      return IDUF

@deprecated
def DATIN(IREAD,DTNORM, efix, ntc, COMS, store,lptfile):
      #return calculate_sample_bins(IREAD,DTNORM, efix, ntc, COMS, store,lptfile)
      
      IDUF=0
      SMALL=1.0E-10
      DSUM=0.0

      # get Average Q value
      COS2TH=2.0*np.cos(COMS["DATA"].theta(IREAD)*pi/180.0)
      QQ=efix+efix-COS2TH*abs(efix)
      COMS["DATA"].QAVRG.set(IREAD, 0.69469*np.sqrt(QQ))
      store.open(53,lptfile)
      store.write(53,' ----------------------------------------------------')
      store.write(53, f' Group {IREAD},  theta = {COMS["DATA"].theta(IREAD)},  Q = {COMS["DATA"].QAVRG(IREAD)}')
      store.write(53,' ----------------------------------------------------')
      store.close(unit=53)
      
      DTNRM=DTNORM(IREAD)
      NDAT=ntc-1

      # laod sample data
      COMS["DATA"].XDAT.copy(COMS["DATA"].xin.output_range(end=COMS["DATA"].NDAT))
      COMS["DATA"].DAT.copy(COMS["DATA"].yin.output_range(end=COMS["DATA"].NDAT))
      # these are already in sample data COM
      for I in get_range(1,COMS["DATA"].NDAT):
       COMS["DATA"].SIG.set(I, COMS["DATA"].ein(I)*DTNRM)
       # only record values if errors are above a tol
       if COMS["DATA"].SIG(I) > SMALL:
        COMS["DATA"].SIG.set(I, pow(COMS["DATA"].SIG(I),2))
        DSUM=DSUM+COMS["DATA"].DAT(I)
       else:
        COMS["DATA"].SIG.set(I, 0.0)
      # if no peak data
      if DSUM < SMALL:
         IDUF=1
      DATIN1(COMS, store, lptfile)
      #
      return IDUF
      


def write_file_info(Nu, ISP, COMS, store, files):
      for n in range(Nu):
       store.open(n, files[n])
       store.write(n,f'{COMS["DATA"].QAVRG(ISP):.3f}, {COMS["SCL"].ASCL:.6f}, {COMS["SCL"].WSCL:.6f}, {COMS["SCL"].BSCL:.6e}, {COMS["SCL"].GSCL:.6e}')
       store.close(unit=n)

@deprecated
def FileInit(Nu, ISP, COMS, store, files):

      for n in range(Nu):
       store.open(n, files[n])
       store.write(n,f'{COMS["DATA"].QAVRG(ISP)}, {COMS["SCL"].ASCL}, {COMS["SCL"].WSCL}, {COMS["SCL"].BSCL}, {COMS["SCL"].GSCL}')
       store.close(unit=n)


def set_sacle_factors(N_QE_peaks,scale,COMS,store, prog,lptfile, o_bgd ):
      """
      This seems to assume that a large error -> its BG measurment
      If its a peak the data should dominate the measurment and the error is small
      Also assume that there are unlikely to be peaks near to the start/end of the data
      """
      # if to scale the values
      if scale <= 1:
         SMALL = 1.0E-10
         y_sum = 0.0
         N_sum = 0
         # get upto 10 "large" values in SIG from first 20 values
         for I in get_range(1,20):
          if COMS["DATA"].SIG(I) >=SMALL:
           N_sum += 1
           y_sum += abs(COMS["DATA"].DAT(I))
          if N_sum >= 10:
              break
         # assume large values dominate and normalise
         BG_scale_1=y_sum/float(N_sum)

         y_sum = 0.0
         N_sum = 0
         # get upto 10 "large" values in SIG from last 20 values
         for I in get_range(1,20):
          if COMS["DATA"].SIG(COMS["DATA"].NDAT-I+1) >= SMALL:
           N_sum += 1
           y_sum += abs(COMS["DATA"].DAT(COMS["DATA"].NDAT-I+1))
          
          if N_sum >= 10:
             break

         # BG scale factors
         COMS["SCL"].BSCL=y_sum/float(N_sum) # average BG value per bin
         if o_bgd == 0:
             COMS["SCL"].BSCL=0.0 # zero background
         elif BG_scale_1<COMS["SCL"].BSCL:
             # use the smaller normalization -> how we account for + or - gradient
            COMS["SCL"].BSCL=BG_scale_1
         COMS["SCL"].BSCL=COMS["SCL"].BSCL/2.0 # this is to account for breaking after 10 recordings and for loop going to 20



         # get the mod largest start/end value
         max_x_value = abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))
         if abs(COMS["DATA"].XDAT(1)) > max_x_value:
            max_x_value=abs(COMS["DATA"].XDAT(1))

         COMS["SCL"].WSCL=max_x_value/3.0 # no idea why 3 -> assume width of peak about 1/6 of the range

         N = 0
         integral_y_data =0.0
         sum_sigma = 0.0 # no idea
         for I in get_range(1,COMS["DATA"].NDAT-1):
          # if not BG
          if COMS["DATA"].SIG(I) >= SMALL:
           N += 1 
           # remove BG and scale by bin width
           integral_y_data += (COMS["DATA"].DAT(I)-COMS["SCL"].BSCL)*(COMS["DATA"].XDAT(I+1)-COMS["DATA"].XDAT(I)) # area of histogram (minus BG)
           sum_sigma += np.sqrt((2.0/COMS["DATA"].SIG(I))) # sigma += sqrt{\frac{2}{sigma(I)} }

         # scale the peak amplitude
         COMS["SCL"].ASCL= float(COMS["DATA"].NDAT)*(integral_y_data/float(N))# integral of the average y value in each bin (when ignoring 0's in the average calc)
         # average error of non-BG data
         sum_sigma = sum_sigma/float(N)
         #       (x range of data)*(sum_sigma/sqrt(N)) 
         sum_sigma=(COMS["DATA"].XDAT(COMS["DATA"].NDAT)-COMS["DATA"].XDAT(1))*sum_sigma/np.sqrt(float(N))

         store.open(53,lptfile)

         # assign scale factors
         # if measurments are less then errors
         if COMS["SCL"].ASCL < sum_sigma:
           store.write(53,' qlm> *** Estimate of Amax is being set to lower bound!')
           store.write(53,f' ( {COMS["SCL"].ASCL:14.7e} --> {SUMSIG:14.7e} )')
           COMS["SCL"].ASCL=sum_sigma
         store.write(53,' ----------------------------------------------------')
         store.close(unit=53)

         COMS["SCL"].ASCL=COMS["SCL"].ASCL/COMS["SCL"].GSCL # divide by average bin width -> dimensions back to the same as y data
         # fit parameters
         COMS["FIT"].FITP.set(1, 1.0) # BG min
         COMS["FIT"].FITP.set(2, 1.0) # BG max
         COMS["FIT"].FITP.set(3, 0.5) # elastic peak amplitude
         # seems to be storing the max values and BG for latter
         for I in get_range(1,2):
           COMS["SCL"].SCLVEC.set(I,1,COMS["SCL"].BSCL)
           COMS["SCL"].SCLVEC.set(I,2,COMS["SCL"].BSCL)
         # store the amplitude
         COMS["SCL"].SCLVEC.set(3,1,COMS["SCL"].ASCL) # elastic peak
         COMS["SCL"].SCLVEC.set(3,2,COMS["SCL"].ASCL)

      COMS["SCL"].SCLVEC.set(4,2,1.0)
      if prog == 'w':
       COMS["SCL"].SCLVEC.set(4,1, COMS["SCL"].ASCL)
       COMS["SCL"].SCLVEC.set(5,2, COMS["SCL"].ASCL)
       COMS["SCL"].SCLVEC.set(6,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
       COMS["SCL"].SCLVEC.set(7,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
      else:
       # fit params for inealstic peaks
       for I in get_range(1,N_QE_peaks):
        COMS["SCL"].SCLVEC.set(3+I,1, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(3+I+I,2, COMS["SCL"].ASCL)
        COMS["SCL"].SCLVEC.set(4+I+I,2, COMS["SCL"].WSCL/COMS["SCL"].GSCL)
      COMS["FIT"].NFEW=0

@deprecated
def PRINIT(NQMAX,IXSCAL,COMS,store, prog,lptfile, o_bgd ):
      """
      This seems to assume that a large error -> its BG measurment
      If its a peak the data should dominate the measurment and the error is small
      """
      #set_sacle_factors
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

         COMS["SCL"].BSCL=SM/float(NSUM) # average value per bin
         # use the smaller normalization -> how we account for + or - gradient
         if BSCL1<COMS["SCL"].BSCL:
            COMS["SCL"].BSCL=BSCL1
         COMS["SCL"].BSCL=COMS["SCL"].BSCL/2.0

         if o_bgd == 0:
             COMS["SCL"].BSCL=0.0 # zero background

         # get the mod largest start/end value
         AXMAX=abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))
         if abs(COMS["DATA"].XDAT(1))>AXMAX:
            AXMAX=abs(COMS["DATA"].XDAT(1))

         COMS["SCL"].WSCL=AXMAX/3.0 # no idea why 3 -> assume width of about 1/6 of the range

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

         # scale the peak amplitude
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
         # fit parameters
         COMS["FIT"].FITP.set(1, 1.0) # BG min
         COMS["FIT"].FITP.set(2, 1.0) # BG max
         COMS["FIT"].FITP.set(3, 0.5) # elastic peak amplitude
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
     

def parameter_error_bars(COVAR, NP):
    # gets the error bars from the diag of the covarience matrix
    SMALL = 1.e-20
    parameter_error = vec(NP)
    for I in get_range(1, NP):
        parameter_error.set(I, sqrt(2.0*abs(COVAR(I,I))+SMALL))
    return parameter_error

@deprecated
def ERRBAR(COVAR,NP):
    # gets the error bars from the diag of the covarience matrix
    SMALL=1.0E-20
    SIGPAR = vec(NP)
    for I in get_range(1,NP):
        SIGPAR.set(I, sqrt(2.0*abs(COVAR(I,I))+SMALL))
    return SIGPAR


def FCTNLG(N):
    A = np.asarray([k+1 for k in range(N)])
    return np.sum(np.log10(A))


def LogLikelihood(COMS, chi2,N_data_points,LOG_HESS_det,N_peaks,NMAX, prog, store, lptfile):
      """   
      log_10 probabilities -> priors
      LOG_HESS_det is the sum of the logs of the diagonal elements of the hessian
      this is an odd way of putting the likelihood (some factors have been dropped such as the ln(10)):
          L = [(4*pi)^N_p/(A*W)^{N_p/2} ]*exp(-(chi^2+D)/2)*PI_i=1,N_P i
          where A and W are the amplitude and width scale factors, D is LOG_HESS_det and N_P is the number of peaks
          the 4*pi is an attempt at normalising the likelihood
          PI signifies the multplication of a series
      """
      CHI_2=chi2*float(N_data_points) # scale to mumber of data points
      LOG_CHI_2=-log10(np.exp(1.))*CHI_2/2.0

      LOG_HESS_det=LOG_HESS_det-float(N_peaks*2)*log10(COMS["SCL"].ASCL*COMS["SCL"].WSCL) # scale the hessian logs correctly

      log_likelihood = LOG_CHI_2 - (0.5*LOG_HESS_det) # the terms from "just the fit"
      # Cost of extra parameters
      log_likelihood -= float(N_peaks)*log10(COMS["SCL"].ASCL*COMS["SCL"].WSCL) # scale factor
      log_likelihood += np.sum(np.log10(np.asarray([k+1 for k in range(N_peaks)]))) # cost of more parameters
      log_likelihood += float(N_peaks)*log10(4.0*pi) #an attempt to normalise prob

      store.open(53, lptfile)
      fit_type = f'{N_peaks} Quasi-elastic lines'

      if prog =='s':
          fit_type = "Stretched exp"
      elif prog == 'w':
          fit_type = "Water"

      store.write(53,f' Log10[Prob({fit_type} |DATA)] = {log_likelihood:11.1f}')
      if N_peaks< NMAX and prog=='l':
         store.write(53,' -------------------------')
      elif N_peaks<1:
         store.write(53,' -------------------------')
            
      store.close(unit=53)
      return log_likelihood, LOG_HESS_det


# log likelyhoods?
def PROBN(COMS, CNORM,NDAT,DETLOG,NFEW,NMAX, prog, store, lptfile):
      # log_10 probabilities -> priors
      CHISQ=CNORM*float(NDAT)
      DETLOG=DETLOG-float(NFEW*2)*log10(COMS["SCL"].ASCL*COMS["SCL"].WSCL)
      CHI2LG=-log10(2.7182818)*CHISQ/2.0

      PROBLG=CHI2LG-(0.5*DETLOG)-(float(NFEW)*log10(COMS["SCL"].ASCL*COMS["SCL"].WSCL))
      PROBLG=PROBLG+FCTNLG(NFEW)
      PROBLG=PROBLG+(float(NFEW)*log10(4.0*pi))

      store.open(53, lptfile)
      if prog=='l':
         store.write(53,f' Log10[Prob({NFEW} Quasi-elastic lines|Data)] = {PROBLG:11.1f}')
         if NFEW< NMAX:
            store.write(53,' -------------------------')
      
      if prog=='s':
          store.write(53,f' Log10[Prob(Stretched exp|Data)] = {PROBLG:11.1f}')
          if NFEW<1:
             store.write(53,' -------------------------')
      
      if prog=='w':
         store.write(53,f' Log10[Prob(Water|Data)] = {PROBLG:11.1f}')
         if NFEW<1:
            store.write(53,' -------------------------')
      
      store.close(unit=53)
      #PRBSV.set(NFEW+1, PROBLG)
      return PROBLG, DETLOG


def PRBOUT(P,NP,NQ):
      #REAL  P(4,m_sp),POUT(4,m_sp)
      POUT = matrix_2(4, m_sp)
      J=NQ
      SM=0.0
      for I in get_range(1,NP):
        SM=SM+pow(10.0, P(I,J)-P(NP,J))
      PLMAX = max([P(3,J), P(4,J)])
      PLNORM=log10(SM)+P(NP,J)
      for I in get_range(1,NP):
        P.set(I,J, P(I,J)-PLMAX)
      
      for I in get_range(1,NP):
        POUT.set(I,J, P(I,J))
      return P, POUT
