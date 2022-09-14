from quasielasticbayes.fortran_python import Vec, Matrix_2D, C_Vec, C_Matrix_2D


class DatCom(object):

    def __init__(self, m_d, m_sp):
        self.XDAT = C_Vec(m_d)
        self.DAT = C_Vec(m_d)
        self.SIG = C_Vec(m_d)
        self.NDAT = 0
        self.xin = C_Vec(m_d)
        self.yin = C_Vec(m_d)
        self.ein = C_Vec(m_d)
        self.theta = C_Vec(m_sp)
        self.QAVRG = C_Vec(m_sp)


class Dintrp(object):
    def __init__(self, m_d):
        self.IPDAT = C_Vec(m_d)  # indicies for mapping to original bins
        self.XPDAT = C_Vec(m_d)  # fractional bin offsets


class ModParams(object):
    def __init__(self):
        self.NBIN = 0
        self.IMIN = 0  # min energy
        self.IMAX = 0  # max energy
        self.RSCL = 0
        self.BNORM = 0


class FFTCom(object):
    def __init__(self, m_d, m_d1, m_d2):
        self.FRES = C_Vec(m_d2)
        self.FWRK = C_Vec(m_d2)
        self.XJ = C_Vec(m_d)
        self.TWOPIK = C_Vec(m_d1)
        self.NFFT = 0


class data_object(object):
    def __init__(self, m_d, m_d1, m_d2):
        self.x_data = C_Vec(m_d)
        self.y_data = C_Vec(m_d)
        self.e_data = C_Vec(m_d)
        self.N = 0
        self.x_bin = C_Vec(m_d)
        self.y_bin = C_Vec(m_d)
        self.e_bin = C_Vec(m_d)
        self.N_bin = 0
        self.phases = C_Vec(m_d1)
        self.FTY = Vec(m_d2, True)
        self.IFTY = Vec(m_d2, True)
        self.N_FT = 0


class FitCom(object):
    def __init__(self, m_d, m_p, m_d1):
        self.FIT = C_Vec(m_d)
        self.RESID = C_Vec(m_d)  # residuals
        self.NFEW = 0  # number inelastic peaks
        self.FITP = C_Vec(m_p)  # paramteters
        self.EXPF = C_Matrix_2D(m_d1, 6)  # exponentials


class GRDCom(object):
    def __init__(self, m_d2, m_d, m_p):
        # components of fit convolved with resolution
        self.DDDPAR = C_Matrix_2D(m_d, m_p)
        # inelastic peaks
        self.FR2PIK = C_Matrix_2D(m_d2, 2)


class ModResidual(object):
    def __init__(self, m_d):
        self.ntr = 0
        self.xres = C_Vec(m_d)
        self.yres = C_Vec(m_d)
        self.eres = C_Vec(m_d)
        self.nrbin = 0
        self.ermin = 0
        self.ermax = 0


class QW1Com(object):
    def __init__(self, m_sp):
        self.QW1 = C_Vec(m_sp)
        self.SIGQW1 = C_Vec(m_sp)
        self.ISPEC = 0


class SCLCom(object):
    # These are all scale factors
    def __init__(self, m_p):
        self.BSCL = 0  # background
        self.ASCL = 0  # amplitude
        self.WSCL = 0  # width
        self.SCLVEC = C_Matrix_2D(m_p, 2)
        self.GSCL = 0  # offset from zero


class STEXP(object):
    def __init__(self):
        self.BETEXP = 1


class WRKCom(object):
    def __init__(self, m_d2):
        self.WORK = Matrix_2D(m_d2, 2)
