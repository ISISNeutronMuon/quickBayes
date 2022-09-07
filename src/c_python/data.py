from quasielasticbayes.fortran_python import vec, matrix_2, c_vec, c_matrix_2


class DatCom(object):

    def __init__(self, m_d, m_sp):
        self.XDAT = c_vec(m_d)
        self.DAT = c_vec(m_d)
        self.SIG = c_vec(m_d)
        self.NDAT = 0
        self.xin = c_vec(m_d)
        self.yin = c_vec(m_d)
        self.ein = c_vec(m_d)
        self.theta = c_vec(m_sp)
        self.QAVRG = c_vec(m_sp)


class Dintrp(object):
    def __init__(self, m_d):
        self.IPDAT = c_vec(m_d)  # indicies for mapping to original bins
        self.XPDAT = c_vec(m_d)  # fractional bin offsets


class ModParams(object):
    def __init__(self):
        self.NBIN = 0
        self.IMIN = 0  # min energy
        self.IMAX = 0  # max energy
        self.RSCL = 0
        self.BNORM = 0


class FFTCom(object):
    def __init__(self, m_d, m_d1, m_d2):
        self.FRES = c_vec(m_d2)
        self.FWRK = c_vec(m_d2)
        self.XJ = c_vec(m_d)
        self.TWOPIK = c_vec(m_d1)
        self.NFFT = 0


class data_object(object):
    def __init__(self, m_d, m_d1, m_d2):
        self.x_data = c_vec(m_d)
        self.y_data = c_vec(m_d)
        self.e_data = c_vec(m_d)
        self.N = 0
        self.x_bin = c_vec(m_d)
        self.y_bin = c_vec(m_d)
        self.e_bin = c_vec(m_d)
        self.N_bin = 0
        self.phases = c_vec(m_d1)
        self.FTY = vec(m_d2, True)
        self.IFTY = vec(m_d2, True)
        self.N_FT = 0


class FitCom(object):
    def __init__(self, m_d, m_p, m_d1):
        self.FIT = c_vec(m_d)
        self.RESID = c_vec(m_d)  # residuals
        self.NFEW = 0  # number inelastic peaks
        self.FITP = c_vec(m_p)  # paramteters
        self.EXPF = c_matrix_2(m_d1, 6)  # exponentials


class GRDCom(object):
    def __init__(self, m_d2, m_d, m_p):
        # components of fit convolved with resolution
        self.DDDPAR = c_matrix_2(m_d, m_p)
        # inelastic peaks
        self.FR2PIK = c_matrix_2(m_d2, 2)


class ModResidual(object):
    def __init__(self, m_d):
        self.ntr = 0
        self.xres = c_vec(m_d)
        self.yres = c_vec(m_d)
        self.eres = c_vec(m_d)
        self.nrbin = 0
        self.ermin = 0
        self.ermax = 0


class QW1Com(object):
    def __init__(self, m_sp):
        self.QW1 = c_vec(m_sp)
        self.SIGQW1 = c_vec(m_sp)
        self.ISPEC = 0


class SCLCom(object):
    # These are all scale factors
    def __init__(self, m_p):
        self.BSCL = 0  # background
        self.ASCL = 0  # amplitude
        self.WSCL = 0  # width
        self.SCLVEC = c_matrix_2(m_p, 2)
        self.GSCL = 0  # offset from zero


class STEXP(object):
    def __init__(self):
        self.BETEXP = 1


class WRKCom(object):
    def __init__(self, m_d2):
        self.WORK = matrix_2(m_d2, 2)
