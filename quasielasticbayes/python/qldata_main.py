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
      GRAD, DPAR = vec(m_p,True), vec(m_p)
      HESS, COVAR = matrix_2(m_p,m_p,True), matrix_2(m_p,m_p)
      DTNORM,XSCALE = vec(m_sp), vec(m_sp)
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
       store.open(1, sfile[:l_fn]+'_test.python.lpt')
       vals =  COMS["FFT"].FRES.output()
       vals2 =  COMS["FFT"].FWRK.output()
       for v in range(2000): 
           store.write(1,str(vals[v])+"  "+ str(vals2[v]))
       store.close(1)



      print("Hi It worked!!!!!!!!!!!!!! #actually doing QL data not res")
      nd_out=10
      store.dump()
      return nd_out,xout.output(),yout.output(),eout.output(),yfit.output(),yprob.output()
