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
      COMS["res_data"] = data_object(m_d, m_d1, m_d2)
      COMS["sample_data"] = data_object(m_d, m_d1, m_d2)
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
      HESS = matrix_2(m_p,m_p)
      PRMSV,SIGSV = matrix_3(7,4,m_sp), matrix_3(7,4,m_sp)
      INDX = vec(m_p)
      LGOOD = BoolVec(m_sp)
      prog='l'

      COMS["FFT"].NFFT=m_d
      #numb = [ngrp, nsp, ntc, Ndat, nbin, IMIN, IMAX, NB, nrbin]
      NSPEC=numb[0] #no. of spectra
      ISP=numb[1] # number of spectra
      COMS["QW1"].ISPEC=ISP
      ntc=numb[2] #no. of points
      COMS["DATA"].NDAT=numb[3]
      COMS["sample_data"].N = int(numb[3])
      COMS["Params"].NBIN=int(numb[4])
      COMS["Params"].IMIN=numb[5]
      COMS["Params"].IMAX=numb[6]
      NB=numb[7]
      COMS["Res"].nrbin=int(numb[8])
      COMS["res_data"].N_bin = int(numb[8])
      #reals = [efix, theta[isp], rscl, bnorm]
      efix=reals[0]
      COMS["DATA"].theta.set(ISP,reals[1])
      COMS["Params"].RSCL=reals[2]                                      
      COMS["Params"].BNORM=reals[3]
      xin =  COMS["DATA"].xin
      yin =  COMS["DATA"].yin
      ein =  COMS["DATA"].ein
      XDAT =  COMS["DATA"].XDAT
      
      # copy sample data
      COMS["DATA"].xin.copy(x_in[0:m_d])
      COMS["DATA"].yin.copy(y_in[0:m_d])
      COMS["DATA"].ein.copy(e_in[0:m_d])
      COMS["DATA"].XDAT.copy(XD_in[0:m_d])

      COMS["sample_data"].x_data.copy(x_in[0:m_d])
      COMS["sample_data"].y_data.copy(y_in[0:m_d])
      COMS["sample_data"].e_data.copy(e_in[0:m_d])

      # copy resolution data
      COMS["res_data"].x_data.copy(X_r[0:NB])
      COMS["res_data"].y_data.copy(Y_r[0:NB])
      COMS["res_data"].e_data.copy(E_r[0:NB])

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

      # set QW values
      if o_w1 ==1:
       COMS["QW1"].QW1.copy(Wy_in[0:NSPEC])
       tmp =  COMS["QW1"].QW1.output_range(start=1,end=NSPEC)
       tmp =0.5*(np.abs(tmp)+0.00001)
       COMS["QW1"].QW1.copy(tmp)

       COMS["QW1"].SIGQW1.copy(We_in[0:NSPEC])
       tmp =  COMS["QW1"].SIGQW1.output_range(start=1,end=NSPEC)
       tmp =0.5*(np.abs(tmp)+0.00001)
       COMS["QW1"].SIGQW1.copy(tmp)

      # report stuff
      if ISP==1: #print info	
       store.open(53,lptfile)
       store.write(53,f' Sample run: {sfile}')
       store.write(53,f' Resolution file: {rfile}')
       store.write(53,f' Energy range: {xin(COMS["Params"].IMIN):8.3f} to {xin(COMS["Params"].IMAX):8.3f} meV')
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
            store.write(n,' Data : ')#+name) # no idea what name is
            store.write(n,title[:l_title])
            store.write(n,user[:l_user])
            store.write(n, f"{NSPEC}     {COMS['DATA'].NDAT}      {xin(COMS['Params'].IMIN):10.3f}     {xin(COMS['Params'].IMAX):10.3f}")
            store.write(n,'-------------------------------------------------')
            store.write(n,rfile)
            store.write(n,n)
            store.write(n,'-------------------------------------------------')
            store.read(n)
            store.close(unit=n)
      IDUF = 0
      XBLR,YBLR=bin_resolution(NB,0,IDUF, COMS, store, lptfile) # rebin + FFT of splined data -> make bins even spaced

      #debug_dump(sfile[:l_fn]+'_test.python.lpt', COMS["FFT"].FRES.output(),  store) # keep this one
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt', flatten(COMS["res_data"].FTY.output()),  store) # keep this one
      bin_offsets(COMS)
      normalize_x_range(COMS) # record fractional original x bins
      # read in sample data and rebin it to even bins
      IDUF = calculate_sample_bins(ISP,DTNORM, efix, ntc, COMS, store, lptfile) 
      # get inverse FFT for resolution data (if IDUF !=0)
      XBLR, YBLR = bin_resolution(NB,ISP,IDUF, COMS, store, lptfile)

      if IDUF!=0:
          LGOOD.set(ISP, False)

      # this might not even do anything
      bin_offsets(COMS) # get offsets for binned sample data
      set_sacle_factors(3,1, COMS, store, prog, lptfile, o_bgd) # rescale the sample data to make it easier to fit. Set up fits

      write_file_info(3,ISP, COMS, store, [fileout1, fileout2, fileout3]) # dump data to file

      DETLOG = 0
      HESS, COVAR, DPAR, DETLOG =refine_param_values(GRAD,HESS, 3+COMS["FIT"].NFEW,DETLOG,INDX,COVAR, COMS, construct_fit_and_chi, prog, o_bgd,o_w1, o_el, store, lptfile)

      ################################
      # chnage the code so it only calculates 0 and 1 elastic lines
      ###############################
      for counter in range(4): # equivalent of less than equal to 3
            if counter > 0: # skip on the first pass -> no peaks to find
               #print("hi")
               HESS, COVAR, DPAR, DETLOG = find_latest_peak(COMS, GRAD,HESS,DPAR, INDX,COVAR, o_w1, prog, o_bgd, o_el, store, lptfile, DETLOG, construct_fit_and_chi)

            NPARMS=4+2*COMS["FIT"].NFEW # 4 default (BG and elastic params) plus 2 per peak
            chi_keep=construct_fit_and_chi(COMS["FIT"].FITP,COMS, o_bgd, o_w1)
            # get copy of fit values
            FITPSV.copy(COMS["FIT"].FITP.output_range(end=NPARMS))

            step_size=0.3
            if COMS["FIT"].NFEW>1: # if more than one peak -> smaller steps
                step_size=step_size/10.0
            IAGAIN=0
            delta_chi_threshold=0.003
            # optimization steps
            for I in get_range(1,200):
                #print("test", I,NPARMS, len(HESS.output()))

                #######################################################################################################
                # come back to this one!!!!!
                HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,step_size, o_bgd, o_w1,o_el, prog)
                #######################################################################################################
                DPAR, fit_params= update_fit_params(COVAR,GRAD,NPARMS,COMS["FIT"].NFEW,COMS["FIT"].FITP,prog,store, lptfile)
                COMS["FIT"].FITP.copy(fit_params)
                CNORM=construct_fit_and_chi(COMS["FIT"].FITP,COMS, o_bgd, o_w1)

                if CNORM<=chi_keep: # a better parameter set
                    #print("a")
                    chi_diff=(chi_keep-CNORM)/CNORM
                    if abs(chi_diff)<=delta_chi_threshold: 
                    #if abs(abs(CHIDIF)-CDIFMN)<=1.e-1: # fudge factor to make sure it does the same as Fortran 
                        #print("waaa")
                        if IAGAIN==0:
                            #print('option1')
                            delta_chi_threshold=0.00005
                            step_size=0.15
                            IAGAIN=1
                        else:
                            #print('go to 3')
                            break
             
                    chi_keep=CNORM
                    FITPSV.copy(COMS["FIT"].FITP.output_range(end=NPARMS))
                else:
                    COMS["FIT"].FITP.copy(FITPSV.output_range(end=NPARMS))
                    step_size=step_size*0.6
                    if step_size<1.0E-10: # step size is too small
                        #print("another go to 3")
                        break
            #######################################################################################################
            # come back to this one!!!!!        
            #if counter ==0:
            #    debug_dump(sfile[:l_fn]+'_test.python2.lpt',HESS.output(),  store) # keep this one

            HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.7, o_bgd, o_w1,o_el, prog)
            #######################################################################################################
            SIGPAR = parameter_error_bars(COVAR,NPARMS)

            tmp_params, tmp_errors = record_fit_results(COMS, SIGPAR,CNORM, store, lptfile)
            PRMSV.copy(tmp_params, 1,COMS["FIT"].NFEW+1,ISP)
            SIGSV.copy(tmp_errors,1,COMS["FIT"].NFEW+1,ISP )
            ##################################################################################################
            OUTPRM(COMS["FIT"].FITP,COVAR,NPARMS,COMS["FIT"].NFEW,CNORM, store, [fileout1, fileout2, fileout3]) # need to translate
        
            HESS, COVAR, DETLOG = REFINE(COMS, GRAD,HESS,NPARMS,DETLOG,INDX,COVAR,0.25, o_bgd, o_w1,o_el, prog)
            #####################################################################################################

            tmp, DETLOG= LogLikelihood(COMS, CNORM,numb[3],DETLOG,COMS["FIT"].NFEW,3, prog, store, lptfile) # use numb encase NDAT in oms has been changed
            PRBSV.set(COMS["FIT"].NFEW+1, ISP, tmp)
            noff=COMS["DATA"].NDAT*COMS["FIT"].NFEW # offset in vector
            for n in get_range(1,COMS["DATA"].NDAT):
                yfit.set(noff+n, COMS["FIT"].FIT(n))
            COMS["FIT"].NFEW+=1
            if o_el==0: # no peak
                COMS["FIT"].FITP.set(3, 0.0)
                FITPSV.set(3, 0.0)
      
      nd_out = COMS["DATA"].NDAT
     
      for n in get_range(1, nd_out):
           xout.set(n, COMS["DATA"].XDAT(n))
           yout.set(n, COMS["DATA"].DAT(n))
           if COMS["DATA"].SIG(n) > 1e-10:
               eout.set(n,np.sqrt(2./COMS["DATA"].SIG(n)))
           else:
               eout.set(n, 0.0)

      PRBSV, POUT = PRBOUT(PRBSV,4,ISP)
      #print(POUT.output())
      for l in get_range(1,4):
           yprob.set(l, POUT(l,ISP))
      debug_dump(sfile[:l_fn]+'_test.python.lpt',HESS.output(),  store) # keep this one
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["Dintrp"].IPDAT.output_range(end=2000),  store) # keep this one
                 
      debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["DATA"].SIG.output(),  store) # keep this one

      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',HESS.output(), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["SCL"].SCLVEC.output_range(1,2,end=2000), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["GRD"].DDDPAR.output_range(1,6,end=2000), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["FIT"].FITP.output(),store)#COMS["FFT"].FWRK.output_range(end=2000), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["FFT"].FWRK.output_range(end=2000), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["Dintrp"].XPDAT.output_range(end=2000), store)
      #debug_dump(sfile[:l_fn]+'_test.python2.lpt',COMS["WORK"].WORK.output_range(1,1,end=2000), store)
      #print("hi", CNORM, DETLOG,PRMSV(1,COMS["FIT"].NFEW+1,ISP),SIGSV(1,COMS["FIT"].NFEW+1,ISP))
      

      #print("Hi It worked!!!!!!!!!!!!!! #actually doing QL data not res")
      store.dump()
      return nd_out,xout.output(),yout.output(),eout.output(),yfit.output(),yprob.output()
