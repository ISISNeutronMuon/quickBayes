# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2022 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +

from  quasielasticbayes.python.fortran_python import vec, matrix_2


class DatCom(object):

    def __init__(self,m_d, m_sp):
        self.XDAT = vec(m_d)
        self.DAT = vec(m_d)
        self.SIG = vec(m_d)
        self.NDAT = 0
        self.xin = vec(m_d)
        self.yin = vec(m_d)
        self.ein = vec(m_d)
        self.theta = vec(m_sp)
        self.QAVRG = vec(m_sp)

class Dintrp(object):
    def __init__(self, m_d):
        self.IPDAT = vec(m_d)
        self.XPDAT = vec(m_d)


class ModParams(object):
    def __init__(self):
        self.NBIN =0
        self.IMIN = 0
        self.IMAX = 0
        self.RSCL = 0
        self.BNORM = 0


class FFTCom(object):
    def __init__(self,m_d, m_d1, m_d2):
        self.FRES = vec(m_d2)
        self.FWRK = vec(m_d2)
        self.XJ = vec(m_d)
        self.TWOPIK = vec(m_d1)
        self.NFFT = 0


class FitCom(object):
    def __init__(self, m_d,m_p,m_d1):
        self.FIT = vec(m_d)
        self.RESID = vec(m_d)
        self.NFEW =0
        self.FITP = vec(m_p)
        self.EXPF = matrix_2(m_d1,6)

                     
class GRDCom(object):
    def __init__(self,m_d2, m_d, m_p):
        self.DDDPAR = matrix_2(m_d,m_p) # covar matrix?
        self.FR2PIK = matrix_2(m_d2,2)

class ModResidual(object):
    def __init__(self,m_d):
        self.ntr = 0
        self.xres = vec(m_d)
        self.yres = vec(m_d)
        self.eres = vec(m_d)
        self.nrbin =0
        self.ermin = 0
        self.ermax = 0

class QW1Com(object):
   def __init__(self, m_sp):
       self.QW1 = vec(m_sp)
       self.SIGQW1 = vec(m_sp)
       self.ISPEC = 0


class SCLCom(object):
    def __init__(self, m_p):
        self.BSCL = 0
        self.ASCL = 0
        self.WSCL = 0
        self.SCLVEC = matrix_2(m_p,2)
        self.GSCL = 0


class STEXP(object):
    def __init__(self):
        self.BETEXP = 1


class WRKCom(object):
   def __init__(self, m_d2):
       self.WORK = matrix_2(m_d2,2)
