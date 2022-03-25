# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from quasielasticbayes.python.fortran_python import *
from quasielasticbayes.python.constants import *
from quasielasticbayes.python.data import *
from quasielasticbayes.python.qldata_subs import *
from quasielasticbayes.python.bayes import *

def QLdata(numb,x_in,y_in,e_in,reals,opft,XD_in,X_r,Y_r,E_r,Wy_in,We_in,sfile,rfile,l_fn):
      COMS = {}
      COMS["FFT"] = FFTCom(m_d, m_d1, m_d2)
      COMS["DATA"] = DatCom(m_d,m_sp)
      COMS["Dintrp"] = Dintrp(m_d)
      COMS["FIT"] = FitCom(m_d,m_p,m_d1)
      COMS["SCL"] = SCLCom(m_p)
      COMS["GRD"] = GRDCom(m_d2, m_d,m_p)
      COMS["QW1"] = QW1Com(m_sp)
      COMS["WORK"] = WRKCom(m_d2)
      COMS["Res"] = ModResidual(m_d)
      COMS["Params"] = ModParams()
      STExp = STEXP()
      store = storage()
      print("PYTHON>>>>>")
      #real x_in(m_d), y_in(m_d), e_in(m_d)
      #intent(in) :: x_in, y_in, e_in                    !sample data
      #real reals(4)
      #intent(in) :: reals                               !real parameters
      #real XD_in(m_d), X_r(m_d), Y_r(m_d), E_r(m_d)
      #intent(in) :: XD_in, X_r, Y_r, E_r                 !sample xrange, res data (blur)
      #real Wy_in(m_sp), We_in(m_sp)
      #intent(in) :: Wy_in, We_in                        !fixed width data
      #integer numb(9)
      #intent(in) :: numb                                !integer parameters
      #integer opft(4)
      #intent(in) :: opft                                 !options parameters
      #integer l_fn
      #intent(in) :: l_fn                                 !length of filenames
      #character*140 sfile, rfile
      #intent(in) :: sfile, rfile                         !sample & res filenames
      nd_out = 0 #number of output points
      xout, yout, eout = vec(m_d), vec(m_d), vec(m_d) #!data values
      yfit = vec(4*m_d) #fit values
      yprob = vec(4) #probability values
    
      XBLR,YBLR = vec(m_d), vec(m_d)
      GRAD = vec(m_p)
      COVAR = matrix_2(m_p,m_p)
      DTNORM,XSCALE = vec(m_sp), vec(m_sp)
      FITPSV = vec(m_p)
      PRBSV,POUT = matrix_2(4,m_sp), matrix_2(4,m_sp)
      PRMSV,SIGSV = matrix_3(7,4,m_sp), matrix_3(7,4,m_sp)
      INDX = vec(m_p)
      LGOOD = BoolVec(m_sp)
      prog='l'

      COMS["FFT"].NFFT=m_d
      #numb = [ngrp, nsp, ntc, Ndat, nbin, IMIN, IMAX, NB, nrbin]
      NSPEC=numb[0] #no. of groups
      ISP=numb[1] #group number
      COMS["QW1"].ISPEC=ISP
      ntc=numb[2] #no. of points
      COMS["DATA"].NDAT=numb[3]
      COMS["Params"].NBIN=int(numb[4])
      COMS["Params"].IMIN=numb[5]
      COMS["Params"].IMAX=numb[6]
      NB=numb[7]
      COMS["Res"].nrbin=int(numb[8])
      #reals = [efix, theta[isp], rscl, bnorm]
      efix=reals[0]
      COMS["DATA"].theta.set(ISP,reals[1])
      COMS["Params"].RSCL=reals[2]                                      
      COMS["Params"].BNORM=reals[3]
      xin =  COMS["DATA"].xin
      yin =  COMS["DATA"].yin
      ein =  COMS["DATA"].ein
      XDAT =  COMS["DATA"].XDAT
      

      COMS["DATA"].xin.copy(x_in[0:m_d])
      COMS["DATA"].yin.copy(y_in[0:m_d])
      COMS["DATA"].ein.copy(e_in[0:m_d])
      COMS["DATA"].XDAT.copy(XD_in[0:m_d])

      COMS["Res"].xres.copy(X_r[0:NB])
      COMS["Res"].yres.copy(Y_r[0:NB])
      COMS["Res"].eres.copy(E_r[0:NB])
      o_el=opft[0]
      o_bgd=opft[1]
      o_w1=opft[2]
      # can speed this up
      DTNORM.fill(1., m_sp) #DTNORM, NSPEC
      XSCALE.fill(1.0, m_sp) #XSCALE, NSPEC
      lptfile=''
      fileout1=''
      fileout2=''
      fileout3=''
      l_lpt=l_fn+8
      lptfile=sfile[:l_fn]+'_QLd.python.lpt'
      l_file=l_fn+8
      fileout1=sfile[:l_fn]+'_QLd.python.ql1'
      fileout2=sfile[:l_fn]+'_QLd.python.ql2'
      fileout3=sfile[:l_fn]+'_QLd.python.ql3'
      l_title=9
      title='<unknown>'
      l_user=9
      user='<unknown>'

      if o_w1 ==1:
       COMS["QW1"].QW1.copy(Wy_in[0:NSPEC])
       tmp =  COMS["QW1"].QW1.output_range(start=1,end=NSPEC)
       tmp =0.5*(np.abs(tmp)+0.00001)
       COMS["QW1"].QW1.copy(tmp)

       COMS["QW1"].SIGQW1.copy(We_in[0:NSPEC])
       tmp =  COMS["QW1"].SIGQW1.output_range(start=1,end=NSPEC)
       tmp =0.5*(np.abs(tmp)+0.00001)
       COMS["QW1"].SIGQW1.copy(tmp)

      if ISP==1: #print info	
       store.open(53,lptfile)
       store.write(53,f' Sample run: {sfile}')
       store.write(53,f' Resolution file: {rfile}')
       store.write(53,f' Energy range: {xin(COMS["Params"].IMIN)} to {xin(COMS["Params"].IMAX)} meV')
       if o_el == 0:
          store.write(53,' Elastic option : NO peak')
       if o_el==1:
          store.write(53,' Elastic option : WITH peak')
       if o_bgd==2:
          store.write(53,' Background option : sloping')
       if o_bgd==1:
          store.write(53,' Background option : flat')
       if o_bgd==0:
          store.write(53,' Background option : zero')

       if o_w1 ==1:
          store.write(53,' Width option : fixed from file ')
       else:
          store.write(53,' Width option : free ')
       
       store.close(unit=53)
       store.open(1,fileout1)
       store.open(2,fileout2)
       store.open(3,fileout3)
       for n in get_range(1,3):
            #store.write(n,' Data : '+name) # no idea what name is
            store.write(n,title[:l_title])
            store.write(n,user[:l_user])
            store.write(n, f"{NSPEC}, {COMS['DATA'].NDAT}, {xin(COMS['Params'].IMIN)}, {xin(COMS['Params'].IMAX)}")
            store.write(n,'-------------------------------------------------')
            store.write(n,rfile)
            store.write(n,n)
            store.write(n,'-------------------------------------------------')
            store.read(n)
            store.close(unit=n)
       IDUF = 0
       XBLR,YBLR=BLRINT(NB,0,IDUF, COMS, store, lptfile) # rebin + FFT of splined data -> improves signal to noise ratio of res file

       #debug_dump(sfile[:l_fn]+'_test.python.lpt', COMS["FFT"].FWRK.output(),  store) # keep this one

       DPINIT(COMS) # subtracts the res from data in time domain
       GDINIT(COMS) # normalised (by GRD) offset of XDAT values -> might be covar matrix
       # seems to filter data -> if large errors assume they dominate the bin
       # sample data into data, XDAT is binned sample data, DAT is binned sample y
       IDUF = DATIN(ISP,DTNORM, efix, ntc, COMS, store, lptfile) 
       # FFT of time binned sample data into FRES
       XBLR, YBLR = BLRINT(NB,ISP,IDUF, COMS, store, lptfile) # back to real space? - print resolution range
       if IDUF!=0:
          LGOOD.set(ISP, False)

       # this might not even do anything
       DPINIT(COMS) # subtraction using sample time domain data
       PRINIT(3,1, COMS, store, prog, lptfile, o_bgd) # seems to find the dominanat data - i.e. not BG and it records the BG
       FileInit(3,ISP, COMS, store, [fileout1, fileout2, fileout3]) # dump data to file
       DETLOG = 0
       HESS, COVAR, DPAR=REFINA(GRAD,3+COMS["FIT"].NFEW,DETLOG,INDX,COVAR, COMS, CCHI, prog, o_bgd,o_w1, o_el, store, lptfile)
       #GOTO 2
       ################################
       # chnage the code so it only calculates 0 and 1 elastic lines
       ###############################
       for counter in range(4): # equivalent of less than equal to 3
            if counter > 0: # skip on the first pass
               print("hi")
               HESS, COVAR, DPAR = SEARCH(COMS, GRAD,HESS,DPAR, INDX,COVAR, o_w1, prog, o_bgd, o_el, store, lptfile, DETLOG, CCHI)
            NPARMS=4+2*COMS["FIT"].NFEW # line 2
            CHIOLD=CCHI(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
            FITPSV.copy(COMS["FIT"].FITP.output_range(end=NPARMS))
            STEPSZ=0.3
            if COMS["FIT"].NFEW>1:
                STEPSZ=STEPSZ/10.0
            IAGAIN=0
            CDIFMN=0.003
            for I in get_range(1,200):
                #print("test", I)

                HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,STEPSZ, o_bgd, o_w1,o_el, prog)

                DPAR = NEWEST(COVAR,GRAD,NPARMS,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)
                CNORM=CCHI(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
                if CNORM<=CHIOLD:
                    #print("a")
                    CHIDIF=(CHIOLD-CNORM)/CNORM
                    if abs(CHIDIF)<=CDIFMN: 
                    #if abs(abs(CHIDIF)-CDIFMN)<=1.e-1: # fudge factor to make sure it does the same as Fortran 
                        print("waaa")
                        if IAGAIN==0:
                            print('option1')
                            CDIFMN=0.00005
                            STEPSZ=0.15
                            IAGAIN=1
                        else:
                            print('go to 3')
                            break
             
                    CHIOLD=CNORM
                    FITPSV.copy(COMS["FIT"].FITP.output_range(end=NPARMS))
                else:
                    COMS["FIT"].FITP.copy(FITPSV.output_range(end=NPARMS))
                    STEPSZ=STEPSZ*0.6
                    if STEPSZ<1.0E-10:
                        print("another go to 3")
                        break
        
            HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.7, o_bgd, o_w1,o_el, prog)
            SIGPAR = ERRBAR(COVAR,NPARMS)

            tmp_p, tmp_s = SEEFIT(COMS, SIGPAR,CNORM, store, lptfile)
            PRMSV.copy(tmp_p, 1,COMS["FIT"].NFEW+1,ISP)
            SIGSV.copy(tmp_s,1,COMS["FIT"].NFEW+1,ISP )
            OUTPRM(COMS["FIT"].FITP,COVAR,NPARMS,COMS["FIT"].NFEW,CNORM, store, [fileout1, fileout2, fileout3])

            HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.25, o_bgd, o_w1,o_el, prog)
            tmp, DETLOG= PROBN(COMS, CNORM,numb[3],DETLOG,COMS["FIT"].NFEW,3, prog, store, lptfile) # use numb encase NDAT in oms has been changed
            PRBSV.set(COMS["FIT"].NFEW+1, ISP, tmp)
            noff=COMS["DATA"].NDAT*COMS["FIT"].NFEW # offset in vector
            for n in get_range(1,COMS["DATA"].NDAT):
                yfit.set(noff+n, COMS["FIT"].FIT(n))
            COMS["FIT"].NFEW+=1
            if o_el==0: # no peak
                COMS["FIT"].FITP.set(3, 0.0)
                FITPSV.set(3, 0.0)
      

       debug_dump(sfile[:l_fn]+'_test.python.lpt',HESS.output(),  store) # keep this one
       #debug_dump(sfile[:l_fn]+'_test.python.lpt',COMS["Dintrp"].IPDAT.output_range(end=2000),  store) # keep this one

       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',HESS.output(), store)
       debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["SCL"].SCLVEC.output_range(1,2,end=2000), store)
       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["GRD"].DDDPAR.output_range(1,6,end=2000), store)
       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["FIT"].FITP.output(),store)#COMS["FFT"].FWRK.output_range(end=2000), store)
       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["FFT"].FWRK.output_range(end=2000), store)
       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["Dintrp"].XPDAT.output_range(end=2000), store)
       #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["WORK"].WORK.output_range(1,1,end=2000), store)
       print("hi", CNORM, DETLOG,PRMSV(1,COMS["FIT"].NFEW+1,ISP),SIGSV(1,COMS["FIT"].NFEW+1,ISP))
      

      print("Hi It worked!!!!!!!!!!!!!! #actually doing QL data not res")
      nd_out=10
      store.dump()
      return nd_out,xout.output(),yout.output(),eout.output(),yfit.output(),yprob.output()
