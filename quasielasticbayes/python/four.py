# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
#from Bayes.fortran_python import *
from  quasielasticbayes.python.constants import *
from math import pi
import numpy as np
import scipy.fft as sc
#from numba import jit
 
#@jit(nopython=True)
#def get_range(start, end, dx=1):
#    return range(start, int(end+1*np.sign(dx)),int(dx)) 

def get_range(start, end, dx=1):
    return range(start, int(end+1*np.sign(dx)),int(dx)) 
def CMPLX(A,B):
    return complex(A,B)

def Re(x):
    y = []
    for k in x:
        y.append(k.real)
    return np.asarray(y)

def flatten(x):
    y =[]
    if any(np.iscomplex(x)):
      for v in range(len(x)):
        y.append(x[v].real)
        y.append(x[v].imag)
    else:
      for v in range(len(x)):
        y.append(x[v].real)
        y.append(0.0)
    return np.asarray(y)

def compress(x):
    y = []
    end = len(x)
    if end%2!=0:
        end -=1
    for k in range(0,end,2):
        y.append(np.complex(x[k], x[k+1]))
    return y
    
#@jit(nopython=True)
def FOUR2_NEG_IFORM (DATA,N,NDIM,ISIGN,IFORM):
    #NDIM is always 1
    N_tot = 1
    N_tot = N_tot*N[0] #N
    N_tot = (N_tot/N[0])*(N[0]/2+1) #N/2+1
    NREM = 1
    JDIM = 1
    IDIM = 1 # 
    NCURR = N[0]
    NCURR=int(N[0]/2)     # N/2
    nn = int(N[0])
    DATA = FOXRL (DATA,nn,NREM,ISIGN,IFORM)  
    N_tot=N_tot/(N[0]/2+1)*N[0] # N[0]
    NPREV=int(N_tot/(N[0]*NREM)) # 1/NREM =1
    DATA=flatten(BOTRV(compress(DATA),NPREV,NCURR,NREM))
    result = COOL2(compress(DATA),NPREV,NCURR,NREM,ISIGN)
    return result
 
#@jit(nopython=True)
def FOXRL (DATA,N,NREM,ISIGN,IFORM):
      TWOPI=2*np.pi*float(ISIGN)
      IP0=2
      IP1=int(IP0*N/2.)
      IP2=int(IP1*NREM) # NREM =1

      J1=IP1+1
      DATA[1]=DATA[J1-1]
      for I2 in get_range(1,IP2,IP1):
        TEMPR=DATA[I2-1]
        DATA[I2-1]=DATA[I2-1]+DATA[I2]
        DATA[I2]=TEMPR-DATA[I2]
      THETA=TWOPI/float(N)
      SINTH=np.sin(THETA/2.)
      ZSTPR=-2.*SINTH*SINTH
      ZSTPI=np.sin(THETA)
      ZR=(1.-ZSTPI)/2.
      ZI=(1.+ZSTPR)/2.
      ZR=1.-ZR
      ZI=-ZI
      I1MIN=IP0+1
      I1MAX=IP0*int(N/4)+1 #2*2+1=5
      for I1 in get_range(I1MIN,I1MAX,IP0):

        for I2 in get_range(I1,IP2,IP1):

            I2CNJ=int(IP0*(N/2+1)-2*I1+I2)#5*(5)-2*1+1
            if I2-I2CNJ >=0:
                if ISIGN*(2*IFORM+1)<0: 
                    DATA[I2]=-DATA[I2]
            else:
                DIFR=DATA[I2-1]-DATA[I2CNJ-1]
                DIFI=DATA[I2]+DATA[I2CNJ]
                TEMPR=DIFR*ZR-DIFI*ZI
                TEMPI=DIFR*ZI+DIFI*ZR
                DATA[I2-1]=DATA[I2-1]-TEMPR
                DATA[I2]=DATA[I2]-TEMPI
                DATA[I2CNJ-1]=DATA[I2CNJ-1]+TEMPR
                DATA[I2CNJ]=DATA[I2CNJ]-TEMPI
                DATA[I2CNJ-1]=2*DATA[I2CNJ-1]
                DATA[I2CNJ]=2*DATA[I2CNJ]
            DATA[I2-1]=2*DATA[I2-1]
            DATA[I2]=2*DATA[I2]
         
        TEMPR=ZR-.5
        ZR=ZSTPR*TEMPR-ZSTPI*ZI+ZR
        ZI=ZSTPR*ZI+ZSTPI*TEMPR+ZI        
       
      # assume N =1
      return DATA
      
#@jit(nopython=True)
def BOTRV(DATA,NPREV,N,NREM):
      IP0=1
      IP1=int(IP0*NPREV)
      IP4=int(IP1*N) 
      IP5=int(IP4*NREM)
      I4REV=1
      I4MAX=IP4
      I1MAX = 1
      IP2 = 1
      for I4 in get_range(1,I4MAX,IP1):
        if I4-I4REV < 0:
            I1MAX=I4+IP1-IP0
            for I1 in get_range(I4,I1MAX,IP0):
                for I5 in get_range(I1,IP5,IP4):
                    I5REV=int(I4REV+I5-I4)
                    TEMP=DATA[I5-1]
                    DATA[I5-1]=DATA[I5REV-1]
                    DATA[I5REV-1]=TEMP
        IP2=int(IP4/2)
        if I4REV-IP2 >0:
            while IP2-IP1 > -1:
                if I4REV-IP2 <=0:
                    break
                I4REV=I4REV-IP2
                IP2=int(IP2/2)
        I4REV=I4REV+IP2
      return DATA
      
#@jit(nopython=True)
def COOL2(DATA,NPREV,N,NREM,ISIGN):
      TWOPI=2*np.pi*float(ISIGN)
      IP0=1
      IP1=int(IP0*NPREV)
      IP4=int(IP1*N)
      IP5=int(IP4*NREM)
      IP2=IP1
      NPART=N # N/2, assume N >> 2 -> line 20
      IP3 = 0
      flag = True
      while NPART -2 >0:
          NPART = int(NPART/4)
          
      flag = True
      if NPART -2 == 0:
        IP3=IP2*2
        I1MAX=IP1
        for I1 in get_range(1,I1MAX,IP0):
            for I5 in get_range(I1,IP5,IP3):
                I3A=I5
                I3B=int(I3A+IP2)
                TMP=DATA[I3B-1]
                #print("p", I1, I5, I3A, I3B, DATA[I3A-1]) # these all match
                DATA[I3B-1]=DATA[I3A-1]-TMP
                DATA[I3A-1]=DATA[I3A-1]+TMP
                flag = True
        # check the IP2= IP3 and if statement! 
        IP2 = IP3
        if IP2-IP4>=0:
            return DATA
      while IP3-IP4<0:
      
            IP3=IP2*4
            THETA=TWOPI/float(IP3/IP1)
            SINTH=np.sin(THETA/2.)
            WSTP=np.cdouble(complex(-2.*SINTH*SINTH,np.sin(THETA)))
            W=1.
            for I2 in get_range(1,IP2,IP1):
                if I2-1 >0:
                    W2=pow(W,2)
                    W3=pow(W,3)
                else:
                    W2 = 1
                    W3 = 1
                I1MAX=I2+IP1-IP0
                for I1 in get_range(I2,I1MAX,IP0):
                    for I5 in get_range(I1,IP5,IP3):
                        I3A=I5
                        I3B=int(I3A+IP2)
                        I3C=int(I3B+IP2)
                        I3D=int(I3C+IP2)
                        #print("p", DATA[I3B-1],I1, I5, "W2", W2,"expect",DATA[I3B-1].imag*W2.real, DATA[I3B-1].real*W2.imag)
                        if I2-1 >0:
                            ##print("hi", I2)
                            DATA[I3B-1]=W2*DATA[I3B-1]
                            DATA[I3C-1]=W*DATA[I3C-1]
                            DATA[I3D-1]=W3*DATA[I3D-1]
                        #print("p2", DATA[I3B-1],I3A, I3B, I3C, I3D)
                        T0=DATA[I3A-1]+DATA[I3B-1]
                        T1=DATA[I3A-1]-DATA[I3B-1]
                        T2=DATA[I3C-1]+DATA[I3D-1]
                        T3=DATA[I3C-1]-DATA[I3D-1]
                        DATA[I3A-1]=T0+T2
                        DATA[I3C-1]=T0-T2
                        TEMP=1j*T3
                        if ISIGN <1:
                            TEMP=-TEMP
                        DATA[I3B-1]=T1+TEMP
                        DATA[I3D-1]=T1-TEMP # here
                W=W*WSTP+W
            IP2=IP3
            #print("p", IP3, IP4, DATA[0], DATA[1], DATA[2], DATA[3], DATA[4])
      return DATA


def FOUR2 (DATA,N,NDIM,ISIGN,IFORM):
      """
      COOLEY-TUKEY FAST FOURIER TRANSFORM IN USASI BASIC FORTRAN.
      MULTI-DIMENSIONAL TRANSFORM, EACH DIMENSION A POWER OF TWO,
      COMPLEX OR REAL DATA.
      TRANSFORM(K1,K2,...) = SUM(DATA(J1,J2,...)*EXP(ISIGN*2*PI*SQRT(-1)
      *((J1-1)*(K1-1)/N(1)+(J2-1)*(K2-1)/N(2)+...))), SUMMED FOR ALL
      J1 AND K1 FROM 1 TO N(1), J2 AND K2 FROM 1 TO N(2),
      ETC. FOR ALL NDIM SUBSCRIPTS.  NDIM MUST BE POSITIVE AND
      EACH N(IDIM) MUST BE A POWER OF TWO.  ISIGN IS +1 OR -1.
      LET NTOT = N(1)*N(2)*...*N(NDIM).  THEN A -1 TRANSFORM
      FOLLOWED BY A +1 ONE (OR VICE VERSA) RETURNS NTOT
      TIMES THE ORIGINAL DATA.  IFORM = 1, 0 OR -1, AS DATA IS
      COMPLEX, REAL OR THE FIRST HALF OF A COMPLEX ARRAY.	TRANSFORM
      VALUES ARE RETURNED TO ARRAY DATA.  THEY ARE COMPLEX, REAL OR
      THE FIRST HALF OF A COMPLEX ARRAY, AS IFORM = 1, -1 OR 0.
      THE TRANSFORM OF A REAL ARRAY (IFORM = 0) DIMENSIONED N(1) BY N(2)
      BY ... WILL BE RETURNED IN THE SAME ARRAY, NOW CONSIDERED TO
      BE COMPLEX OF DIMENSIONS N(1)/2+1 BY N(2) BY ....  NOTE THAT IF
      IFORM = 0 OR -1, N(1) MUST BE EVEN, AND ENOUGH ROOM MUST BE
      RESERVED.  THE MISSING VALUES MAY BE OBTAINED BY COMPLEX CONJUGA-
      TION.	THE REVERSE TRANSFORMATION, OF A HALF COMPLEX ARRAY DIMEN-
      SIONED N(1)/2+1 BY N(2) BY ..., IS ACCOMPLISHED BY SETTING IFORM
      TO -1.  IN THE N ARRAY, N(1) MUST BE THE TRUE N(1), NOT N(1)/2+1.
      THE TRANSFORM WILL BE REAL AND RETURNED TO THE INPUT ARRAY.
      RUNNING TIME IS PROPORTIONAL TO NTOT*LOG2(NTOT), RATHER THAN
      THE NAIVE NTOT**2.  FURTHERMORE, LESS ERROR IS BUILT UP.
      WRITTEN BY NORMAN BRENNER OF MIT LINCOLN LABORATORY, JUNE 1968. 
      SEE-- IEEE AUDIO TRANSACTIONS (JUNE 1967), SPECIAL ISSUE ON FFT.
      """
      tmp = DATA.output()
      if any(np.iscomplex(tmp)):
        tmp=flatten(tmp)
      result = []
      if IFORM==-1:
          # N < len(DATA) since we dont actually use all of the data
          result = FOUR2_NEG_IFORM(tmp,[N],NDIM,ISIGN,IFORM) # This seems to take data of the form exp(-t*gamma)*cos(omega*t) and pick out the omega freq
      elif ISIGN ==1 and IFORM ==1:
          result = np.conj(sc.fft(tmp))
      elif ISIGN ==-1  and IFORM ==1:
          result = sc.fft(tmp)
      elif ISIGN ==1 and IFORM ==0:
          result = np.conj(sc.fft(tmp[0:N]))
          result = result[0:int(N/2)+1]
      elif ISIGN ==-1  and IFORM ==0:
          result = sc.fft(tmp[0:N])
          result = result[0:int(N/2)+1]
      return np.asarray(result)