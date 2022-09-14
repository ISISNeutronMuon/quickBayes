from quasielasticbayes.fortran_python import (
    update_vec,
    round_sig,
    get_range,
    NINT,
    Vec,
    Matrix_2D,
    find_index,
    deprecated)

from quasielasticbayes.constants import m_d, m_d2

from quasielasticbayes.four_python import (
    flatten, compress, FOUR2, FOUR3, FOUR2_raw)

from quasielasticbayes.bayes import (
    construct_gradients,
    complex_shift,
    VMLTRC,
    VMLTIC,
    bin_shift,
    VRDOTR,
    HESS0,
    make_hessian,
    INVERT,
    refine_param_values)

from math import sqrt
import numpy as np
from scipy.interpolate import interp1d


def SPLINE(X, Y, N, YP1, YPN, Y2):
    """
    X,Y are vals
    YP1 is the derivative at start
    YPN is the '' at end
    Y2 is the list of second derivatives
    """
    return interp1d(X.output(), Y.output(), kind='cubic')


def SPLINT(X, func):
    return func(X)


def rm_BG(x_bin, y_bin, N_bin, YMAX, LST, COMS, store, lptfile):
    XDMIN = COMS["DATA"].XDAT(1)
    XDMAX = COMS["DATA"].XDAT(COMS["DATA"].NDAT)
    Y0 = YMAX / 10.0
    # these two check seem to be finiding the first/last y value
    # greater than some ref value
    # that is also in the x range

    def check_data_from_below(XB, YB, Y0, XDMIN, II):
        return YB(II) >= Y0 and XB(II) > XDMIN

    II = find_index((x_bin, y_bin, Y0, XDMIN), 1, N_bin, check_data_from_below)
    XMIN = x_bin(II)

    def check_data_from_above(XB, YB, Y0, XDMAX, II):
        return YB(II) >= Y0 and XB(II) < XDMAX

    II = find_index((x_bin, y_bin, Y0, XDMAX), N_bin,
                    1, check_data_from_above, step=-1)
    XMAX = x_bin(II)

    # get x range of interesting values
    data_window = XMAX - XMIN
    bin_width = data_window / 20.0

    AXMAX = abs(COMS["DATA"].XDAT(1))
    # assume data is approximatly symmetrical
    if abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT)) > AXMAX:
        AXMAX = abs(COMS["DATA"].XDAT(COMS["DATA"].NDAT))
    XNDMAX = 500.0
    if COMS["DATA"].NDAT > int(XNDMAX):
        XNDMAX = float(COMS["DATA"].NDAT)

    mean_bin_width = 2.0 * AXMAX / XNDMAX  # average bin width
    if mean_bin_width > bin_width:  # use the largest bin width
        bin_width = mean_bin_width

    XNGD = (2.0 * AXMAX) / bin_width  # average bin size
    NGD = NINT(np.log(XNGD - 1.0) / np.log(2.0)) + 1
    NGD = pow(2, NGD)  # number of data points

    if NGD > m_d:
        store.open(53, lptfile)
        store.write(53, ' ERROR in XGINIT : too many points')
        store.close(unit=53)
        return
    COMS["FFT"].NFFT = NGD
    COMS["res_data"].N_FT = NGD // 2

    # Set the x values for evenly spaced bins
    COMS["FFT"].XJ.set(1, - bin_width * float(COMS["FFT"].NFFT / 2))
    for j in get_range(2, COMS["FFT"].NFFT):
        COMS["FFT"].XJ.set(j, COMS["FFT"].XJ(j - 1) + bin_width)
    # get the first and last bin
    XMIN = XMIN - 5.0 * data_window
    XMAX = XMAX + 5.0 * data_window
    if XMIN < x_bin(1):
        XMIN = x_bin(1)
    if XMAX > x_bin(N_bin):
        XMAX = x_bin(N_bin)
    if XMIN < XDMIN:
        XMIN = XDMIN
    if XMAX > XDMAX:
        XMAX = XDMAX
    if LST:
        store.open(53, lptfile)

        store.write(53, f' Resolution Range: {XMIN:7.1f} to {XMAX:7.1f} ueV')
        store.close(unit=53)

    # get x range -> via indices
    def check_data_x_min(XB, XMIN, II):
        return XB(II) >= XMIN
    II = find_index((x_bin, XMIN), 1, N_bin, check_data_x_min)
    IMIN = II

    def check_data_x_max(XB, XMAX, II):
        return XB(II) <= XMAX

    II = find_index((x_bin, XMAX), N_bin, 1, check_data_x_max, step=-1)
    IMAX = II

    BG1 = 0.0
    BG2 = 0.0
    # get mean value for 5 bins (that are in range) closest to min/max
    for II in get_range(1, 5):
        BG1 = BG1 + y_bin(IMIN + II - 1)
        BG2 = BG2 + y_bin(IMAX - II + 1)
    BG1 = BG1 / 5.0
    BG2 = BG2 / 5.0

    # the gradiant of the BG
    # no idea where the 4 comes from
    DB = (BG2 - BG1) / float(max(IMAX - IMIN - 4, 1))
    BG = BG1
    # Remove BG value from data
    for II in get_range(IMIN, IMAX):
        y_bin.set(II, y_bin(II) - BG)
        BG = BG + DB

    return XMIN, XMAX, y_bin


@deprecated
def XGINIT(XB, YB, NB, YMAX, LST, COMS, store, lptfile):
    return rm_BG(XB, YB, NB, YMAX, LST, COMS, store, lptfile)


"""
***<set up blur function>**********************************************
"""


def rebin(x_in, y_in, e_in, N_new_bins, N_merged_bins):
    XB = Vec(N_new_bins)
    YB = Vec(N_new_bins)
    SMALL = 1.0E-20
    BNORM = 1.0 / float(N_merged_bins)
    N = 0
    for II in get_range(1, N_new_bins, N_merged_bins):  # new binning
        N += 1
        x_value = 0.0
        y_value = 0.0
        K = 0
        for J in get_range(0, N_merged_bins - 1):  # loop over bins in new bin
            IJ = II + J
            if IJ <= N_new_bins:

                x_value += x_in(IJ)
                if e_in(IJ) > SMALL:  # only include non-zero errors
                    K = K + 1
                    y_value += y_in(IJ)

            XB.set(N, BNORM * x_value)
            YB.set(N, 0.0)
            if K > 0:
                YB.set(II, BNORM * y_value)  # normalise data
    return XB, YB


@deprecated
def BINBLR(WX, WY, WE, NB, NBIN):
    rebin(WX, WY, WE, NB, NBIN)


def bin_resolution(N_bin, IREAD, IDUF, COMS, store, lptfile):
    LSTART = True
    if IREAD == 0:
        LSTART = False

    SMALL = 1.0E-20
    COMS["FFT"].NFFT = m_d
    COMS["res_data"].N_FT = m_d // 2
    # copy resolution data
    xr = Vec(m_d)
    yr = Vec(m_d)
    er = Vec(m_d)
    xr.copy(COMS["Res"].xres.output_range(1, N_bin))
    yr.copy(COMS["Res"].yres.output_range(1, N_bin))
    er.copy(COMS["Res"].eres.output_range(1, N_bin))

    # rebin the resolution data as it is not on an evenly spaced grid
    x_bin, y_bin = rebin(xr, yr, er, N_bin, COMS["Res"].nrbin)

    # sample data x range
    XDMIN = COMS["DATA"].XDAT(1)
    XDMAX = COMS["DATA"].XDAT(COMS["DATA"].NDAT)

    # get indicies for valid x range
    first_index = np.argmax(x_bin.output() >= XDMIN) + 1
    last_index = np.argmin(x_bin.output() <= XDMAX) + 1
    # get total and max Y binned values within valid x range
    YSUM = np.sum(y_bin.output_range(first_index, last_index))
    YMAX = np.max(y_bin.output_range(first_index, last_index))
    if YSUM < SMALL:
        store.open(53, lptfile)
        store.write(53, ' Ysum is too small')
        IDUF = 1
        store.close(unit=53)
        return x_bin, y_bin, IDUF

    XBMIN, XBMAX, y_bin = rm_BG(
        x_bin, y_bin, N_bin, YMAX, LSTART, COMS, store, lptfile)

    # populate FRES with spline of evenly spaced binned data -> data to FFT
    # later
    DER2 = Vec(m_d)
    data = Vec(COMS["FFT"].NFFT)
    func = SPLINE(x_bin, y_bin, N_bin, 0.0, 0.0, DER2)
    # factor of 1/2 "hidden" in N_FT compared to NFFT
    TWOPIN = np.pi / float(COMS["res_data"].N_FT)

    # clean up old data
    COMS["FFT"].FRES.fill(0.0, COMS["FFT"].NFFT)
    COMS["res_data"].FTY.fill(0.0, COMS["res_data"].N_FT)

    # spline the data ontop the sample data bins
    XX = 0.0
    bin_width = COMS["FFT"].XJ(2) - COMS["FFT"].XJ(1)
    data = Vec(COMS["FFT"].NFFT)
    data.set(1, SPLINT(XX, func))

    # use symmetry to have the range we need to loop
    for II in get_range(1, int(COMS["FFT"].NFFT / 2)):
        XX += bin_width
        if XX < XBMAX:
            data.set(II + 1, SPLINT(XX, func))
        if -XX > XBMIN:

            data.set(COMS["FFT"].NFFT + 1 - II, SPLINT(-XX, func))
        # set phases
        COMS["FFT"].TWOPIK.set(II, TWOPIN * float(II - 1)
                               )  # looks to be the phase
        COMS["res_data"].phases.set(II, TWOPIN * float(II - 1))
    COMS["FFT"].TWOPIK.set(int(COMS["FFT"].NFFT / 2) + 1,
                           TWOPIN * float(COMS["FFT"].NFFT / 2))
    COMS["res_data"].phases.set(
        int(COMS["res_data"].N_FT) + 1, TWOPIN * float(COMS["res_data"].N_FT))

    # average y value in each bin -> integral of all data is 1
    SUM = np.sum(data.output())
    BNORM = 1. / (SUM * float(COMS["FFT"].NFFT))
    data.copy(data.output() * BNORM)  # scale the splined data
    out = FOUR2(data, COMS["FFT"].NFFT, 1, 1, 0)
    COMS["FFT"].FRES.copy(flatten(out))
    COMS["res_data"].FTY.copy(out)

    # some rotations?
    for II in get_range(3, COMS["FFT"].NFFT, 4):
        COMS["FFT"].FRES.set(II, -COMS["FFT"].FRES(II))
        COMS["FFT"].FRES.set(II + 1, -COMS["FFT"].FRES(II + 1))

    for k in get_range(2, COMS["res_data"].N_FT, 2):
        COMS["res_data"].FTY.set(k, -COMS["res_data"].FTY(k))

    # IFT of resolution
    if not LSTART:
        tmp = COMS["FFT"].FRES.output_range(end=COMS["FFT"].NFFT + 2)
        COMS["FFT"].FWRK.copy(tmp)
        out = FOUR2(COMS["FFT"].FWRK, COMS["FFT"].NFFT, 1, -1, -1)
        COMS["FFT"].FWRK.copy(flatten(out[0:m_d2]))

        COMS["res_data"].IFTY.copy(
            COMS["res_data"].FTY.output_range(
                end=COMS["res_data"].N_FT))
        # ##### check this ######
        out2 = FOUR3(COMS["res_data"].IFTY, COMS["res_data"].N_FT, 1, -1, -1)
        COMS["res_data"].IFTY.copy(out2[0:m_d2])

    LSTART = True
    return x_bin, y_bin, IDUF


@deprecated
def BLRINT(NB, IREAD, IDUF, COMS, store, lptfile):
    return bin_resolution(NB, IREAD, IDUF, COMS, store, lptfile)


def construct_fit_and_chi(fit_params, COMS, o_bgd, o_w1):
    CHI = 0.0
    min_BG = COMS["SCL"].BSCL * fit_params(1)  # BG 1
    max_BG = COMS["SCL"].BSCL * fit_params(2)  # BG 2
    elastic_amplitude = COMS["SCL"].ASCL * \
        fit_params(3)  # elastic peak amplitude
    av_bin_width = fit_params(4)  # average bin width
    NFT2 = 1 + COMS["FFT"].NFFT // 2

    # get resolution values
    FT_res = COMS["res_data"].FTY.output_range(end=NFT2)
    two_pi_k = COMS["res_data"].phases.output_range(end=NFT2)

    RKEXP, RKEXP2 = complex_shift(FT_res, av_bin_width, two_pi_k)

    COMS["GRD"].FR2PIK.copy(flatten(RKEXP))  # lets make these complex later
    COMS["GRD"].FR2PIK.copy(flatten(RKEXP2), 1, 2)

    # uniform value equal to estimate for elastic peak
    fit_values = np.full(NFT2 + 1, elastic_amplitude)

    # estimate for x bin width?
    XNSCL = -np.log(1.0E-7) / (two_pi_k[1] - two_pi_k[0])

    for J in get_range(1, COMS["FIT"].NFEW):
        amplitude_J = COMS["SCL"].ASCL * fit_params(3 + J + J)
        sigma_J = COMS["SCL"].WSCL * fit_params(4 + J + J) / COMS["SCL"].GSCL
        # this might need a plus 1 to N later
        COMS["FIT"].EXPF.fill(0.0, NFT2, 1, J)
        # estimate for range that the exp will alter the value
        NXFT2 = 1 + NINT(XNSCL / (np.abs(sigma_J) + 1.0E-10))
        if NXFT2 > NFT2:
            NXFT2 = NFT2

        EXP_IJ = np.exp(-COMS["FFT"].TWOPIK.output_range(end=NXFT2) * sigma_J)
        COMS["FIT"].EXPF.copy(EXP_IJ, 1, J)
        fit_values[0:NXFT2 + 1] += amplitude_J * EXP_IJ

    # this relies on having the fit values
    COMS["WORK"].WORK.copy(fit_values, 1, 1)

    # multiply fit by FT(resolution)*(other stuff)
    fit_values = fit_values * RKEXP

    COMS["sample_data"].FTY.copy(fit_values)

    # not sure why we faltten the data
    fit_in_original_domain = flatten(
        FOUR2(COMS["sample_data"].FTY, COMS["FFT"].NFFT, 1, -1, -1))
    fit_in_original_domain = bin_shift(
        fit_in_original_domain,
        COMS,
        1)  # shift back to original binning
    COMS["FIT"].FIT.copy(fit_in_original_domain)

    X1 = COMS["sample_data"].x_data(1)
    BG_grad = 0.0
    if o_bgd == 2:
        BG_grad = (max_BG - min_BG) / \
            (COMS["sample_data"].x_data(COMS["sample_data"].N) - X1)

    xdat = COMS["sample_data"].x_data.output_range(end=COMS["DATA"].NDAT)
    data = COMS["DATA"].DAT.output_range(end=COMS["DATA"].NDAT)
    weights = COMS["DATA"].SIG.output_range(end=COMS["DATA"].NDAT)

    fit_in_original_domain += min_BG + BG_grad * (xdat - X1)  # add BG back in

    diff = fit_in_original_domain - data
    resid = diff * weights  # residuals
    CHI = 0.0

    for K in range(COMS["DATA"].NDAT):
        resid[K] = round_sig(round_sig(diff[K]) * round_sig(weights[K]))
        CHI += diff[K] * resid[K]

    COMS["FIT"].FIT.copy(fit_in_original_domain)
    COMS["FIT"].RESID.copy(resid)

    if o_w1 == 1 and COMS["FIT"].NFEW >= 1:
        RESW1D = (COMS["SCL"].WSCL * fit_params(6) - COMS["QW1"].QW1(
            COMS["QW1"].ISPEC)) / COMS["QW1"].SIGQW1(COMS["QW1"].ISPEC)
        CHI = CHI + 2.0 * pow(RESW1D, 2)
    return CHI / (2.0 * float(COMS["DATA"].NDAT))


@deprecated
def CCHI(V, COMS, o_bgd, o_w1):
    return construct_fit_and_chi(V, COMS, o_bgd, o_w1)


def refine_matrices(COMS, GRAD, HESS, NP,
                    DETLOG, INDX, COVAR, STEPSZ,
                    o_bgd, o_w1, o_el, prog):

    NFT2 = int(COMS["FFT"].NFFT / 2) + 1
    NDAT = COMS["DATA"].NDAT
    NFFT = COMS["FFT"].NFFT
    _ = construct_fit_and_chi(COMS["FIT"].FITP, COMS, o_bgd, o_w1)
    HESS = Matrix_2D(NP, NP)

    # update freq copy
    freq_copy = COMS["FFT"].FWRK.output()

    tmp = VMLTRC(
        COMS["WORK"].WORK.output_range(
            1, 1, end=NFT2 + 1),
        compress(
            COMS["GRD"].FR2PIK.output_range(
                1, 2, end=2 * NFT2 + 2)))  # multiple work by complex vectors

    freq_copy = update_vec(flatten(tmp), freq_copy)
    tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
    IFT_copy = update_vec(flatten(tmp), freq_copy)  # IFT = inverse FT
    # map inverse transform back onto "real" space
    COMS["GRD"].DDDPAR.copy(bin_shift(IFT_copy, COMS), 1, 4)

    # update work copy
    work_copy = COMS["WORK"].WORK.output()
    tmp = VMLTRC(COMS["FFT"].TWOPIK.output_range(
        end=NFT2 + 1), compress(IFT_copy[:2 * (NFT2 + 2)]))

    work_copy = update_vec(flatten(tmp), work_copy)  # COMS
    tmp = compress(work_copy[:2 * (NFT2 + 1)])
    tmp = flatten(VMLTIC(tmp))
    work_copy = update_vec(tmp, work_copy)
    tmp = FOUR2_raw(work_copy, NFFT, 1, -1, -1)
    work_copy = update_vec(flatten(tmp), work_copy)
    COMS["WORK"].WORK.copy(work_copy)

    IFT_copy = update_vec(bin_shift(work_copy, COMS), IFT_copy)
    HESS.set(4, 4, VRDOTR(COMS["FIT"].RESID.output_range(
        end=NDAT), IFT_copy[:NDAT + 1], NDAT - 1))
    freq_copy = update_vec(
        COMS["GRD"].FR2PIK.output_range(
            1, 1, end=NFFT + 2), freq_copy)

    tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
    IFT_copy = update_vec(flatten(tmp), IFT_copy)
    COMS["GRD"].DDDPAR.copy(bin_shift(IFT_copy, COMS), 1, 3)
    freq_copy = update_vec(
        COMS["GRD"].FR2PIK.output_range(
            1, 2, end=NFFT + 2), freq_copy)
    tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
    IFT_copy = update_vec(flatten(tmp), IFT_copy)

    work_copy = update_vec(bin_shift(IFT_copy, COMS), work_copy)
    tmp = VRDOTR(COMS["FIT"].RESID.output_range(
        end=NDAT), work_copy[:NDAT + 1], NDAT - 1)
    HESS.set(3, 4, tmp)
    HESS.set(4, 3, tmp)

    for II in get_range(1, COMS["FIT"].NFEW):
        J = 3 + II + II
        AJ = COMS["FIT"].FITP(J) * COMS["SCL"].ASCL

        tmp = VMLTRC(COMS["FIT"].EXPF.output_range(1, II, end=NFT2 + 2),
                     compress(COMS["GRD"].FR2PIK.output_range(1, 1,
                                                              2 * (2 + NFT2))))

        tmp_freq_copy = update_vec(
            flatten(tmp), work_copy)  # it needs this to work

        # construct freq * exp stuff
        freq_copy = update_vec(flatten(tmp), freq_copy)

        tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
        IFT_copy = update_vec(flatten(tmp), IFT_copy)

        # transform back onto original bins
        COMS["GRD"].DDDPAR.copy(bin_shift(IFT_copy, COMS), 1, J)

        # rotate freq values
        tmp = VMLTRC(COMS["FFT"].TWOPIK.output_range(end=NFT2),
                     compress(tmp_freq_copy[:2 * (NFT2 + 1)]))
        freq_copy = update_vec(flatten(tmp), freq_copy)

        tmp_freq_copy = update_vec(freq_copy[:NFFT + 3], tmp_freq_copy)

        tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
        IFT_copy = update_vec(flatten(tmp), IFT_copy)
        COMS["GRD"].DDDPAR.copy(bin_shift(freq_copy, COMS), 1, J + 1)

        tmp, HESS = HESS0(
            HESS, COMS["FIT"].RESID.output_range(
                end=NDAT), COMS["GRD"].DDDPAR.output_range(
                1, J + 1, NDAT + 1), AJ, J, NDAT)

        COMS["GRD"].DDDPAR.copy(tmp, 1, J + 1)

        tmp = compress(tmp_freq_copy[:2 * (NFT2 + 2)])
        tmp_freq_copy = update_vec(
            flatten(
                VMLTIC(tmp)),
            tmp_freq_copy)  # freq*i

        freq_copy = update_vec(tmp_freq_copy[:NFFT + 2], freq_copy)

        tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
        IFT_copy = update_vec(flatten(tmp), IFT_copy)

        work_copy_col_2 = bin_shift(IFT_copy, COMS)

        tmp = VRDOTR(
            COMS["FIT"].RESID.output_range(
                end=COMS["DATA"].NDAT - 1),
            work_copy_col_2, NDAT - 1)

        HESS.set(4, J, tmp)
        HESS.set(J, 4, tmp)

        tmp = VMLTRC(COMS["FFT"].TWOPIK.output_range(
            end=NFT2), compress(work_copy[:2 * (NFT2 + 1)]))
        freq_copy = update_vec(flatten(tmp), freq_copy)
        COMS["WORK"].WORK.copy(freq_copy[:NFFT + 2])

        tmp = FOUR2_raw(freq_copy, NFFT, 1, -1, -1)
        IFT_copy = update_vec(flatten(tmp), IFT_copy)

        work_copy_col_2 = update_vec(
            bin_shift(IFT_copy, COMS), work_copy_col_2)

        # residuals mapped onto origianl bins
        SM = VRDOTR(COMS["FIT"].RESID.output_range(
            end=NDAT - 1), work_copy_col_2[:NDAT], NDAT - 1)
        HESS.set(4, J + 1, -AJ * SM)
        HESS.set(J + 1, 4, -AJ * SM)

        # this bit is odd, the dimensions don't follow but it works (assuming
        # FOUR2_raw is an IFT)
        tmp = compress(work_copy[:2 * (NFT2)])
        work_copy = update_vec(flatten(VMLTIC(tmp)), work_copy)
        tmp = FOUR2_raw(work_copy, NFFT, 1, -1, -1)
        work_copy = update_vec(flatten(tmp), work_copy)
        freq_copy = update_vec(bin_shift(work_copy, COMS), freq_copy)
        COMS["FFT"].FWRK.copy(bin_shift(work_copy, COMS), 1)

        SM = VRDOTR(COMS["FIT"].RESID.output_range(
            end=NDAT), freq_copy[:NDAT + 1], NDAT)

        HESS.set(J + 1, J + 1, -AJ * SM)

    GRAD.copy(
              construct_gradients(COMS["FIT"].RESID, NDAT, NP,
                                  COMS["SCL"].SCLVEC, COMS, col=2))

    HESS = make_hessian(NP, COMS["SCL"].SCLVEC.output_col(2),
                        STEPSZ, COMS["FIT"].NFEW, prog, COMS,
                        o_el, HESS)

    if o_w1 == 1 and NP > 6:
        DIF = COMS["SCL"].WSCL * \
            COMS["FIT"].FITP(6) - COMS["QW1"].QW1(COMS["QW1"].ISPEC)
        SIG2 = 2.0 / pow(COMS["QW"].SIGQW1(COMS["QW"].ISPEC), 2)
        GRAD.set(6, GRAD(6) + SIG2 * DIF * COMS["SCL"].WSCL)
        HESS.set(6, 6, HESS(6, 6) + SIG2 * pow(COMS["SCL"].WSCL, 2))
    covar_default = 1.
    HESS, COVAR, DETLOG = INVERT(NP, INDX, covar_default, HESS)
    return HESS, COVAR, DETLOG


def REFINE(COMS, GRAD, HESS,
           NP, DETLOG, INDX,
           COVAR, STEPSZ,
           o_bgd, o_w1, o_el, prog):

    return refine_matrices(COMS, GRAD, HESS,
                           NP, DETLOG, INDX,
                           COVAR, STEPSZ,
                           o_bgd, o_w1, o_el,
                           prog)


def record_fit_results(COMS, SIGPAR, Chi2, store, lptfile):
    store.open(53, lptfile)
    scale_error = sqrt(Chi2)
    parameters = []
    parameter_errors = []

    store.write(
        53,
        f' Best-fit assuming no. of quasi-elastic lines = { COMS["FIT"].NFEW}')
    store.write(53, f' >>> Normalised Chi-squared = {Chi2:11.4f}')
    store.write(
        53,
        f' Background(Xmin) = {COMS["FIT"].FITP(1)*COMS["SCL"].BSCL:12.3e} +- {SIGPAR(1)*COMS["SCL"].BSCL*scale_error:12.3e}')  # noqa E501
    store.write(
        53,
        f' Background(Xmax) = {COMS["FIT"].FITP(2)*COMS["SCL"].BSCL:12.3e}  +- {SIGPAR(2)*COMS["SCL"].BSCL*scale_error:12.3e}')  # noqa E501
    store.write(
        53,
        f' Zero offset       = {COMS["FIT"].FITP(4)*COMS["SCL"].GSCL*1000.:12.2f}  +- {SIGPAR(4)*COMS["SCL"].GSCL*scale_error*1000.:10.2f}  ueV')  # noqa E501
    store.write(53, ' Elastic line')

    parameters.append(COMS["FIT"].FITP(3) * COMS["SCL"].ASCL)
    parameter_errors.append(SIGPAR(3) * COMS["SCL"].ASCL * scale_error)
    store.write(
        53,
        f'Amplitude  =   {parameters[-1]:13.4e}  +- {parameter_errors[-1]:11.3e} ')  # noqa E501

    for II in get_range(1, COMS["FIT"].NFEW):
        J = 4 + II + II
        store.write(53, f' Quasi-elastic line {II}')

        parameters.append(2.0 * COMS["FIT"].FITP(J) * COMS["SCL"].WSCL)
        parameter_errors.append(
            2.0 *
            SIGPAR(J) *
            COMS["SCL"].WSCL *
            scale_error)
        # convert from nev to uev
        store.write(
            53,
            f' FWHM        =   {1000*parameters[-1]:13.2f}  +- {1000*parameter_errors[-1]:11.2f} ueV')  # noqa E501

        parameters.append(COMS["FIT"].FITP(J - 1) * COMS["SCL"].ASCL)
        parameter_errors.append(SIGPAR(J - 1) * COMS["SCL"].ASCL * scale_error)
        store.write(
            53,
            f' Amplitude  =   {parameters[-1]:13.4e}  +-  {parameter_errors[-1]:11.3e}')  # noqa E501

    store.close(53)
    return np.asarray(parameters), np.asarray(parameter_errors)


@deprecated
def SEEFIT(COMS, SIGPAR, CNORM, store, lptfile):
    return record_fit_results(COMS, SIGPAR, CNORM, store, lptfile)


def OUTPRM(P, C, NP, NFEW, CNORM, store, files):
    if NFEW < 1 or NFEW > len(files):
        return
    for k in range(len(files)):
        store.open(k + 1, files[k])
    store.write(
        NFEW,
        f'{P(3):13.4e}   {P(1):13.4e}   {P(2):13.4e}   {P(4):13.4e}')

    for II in get_range(5, NP, 2):
        store.write(NFEW, f'{P(II):13.4e}   {P(II+1):13.4e}')

    CSCALE = 2.0 * CNORM
    for J in get_range(1, NP):
        for II in get_range(1, NP):
            C.set(II, J, CSCALE * C(II, J))
    store.write(NFEW, f'{C(3,3):13.4e}')
    store.write(NFEW, f'{C(3,5):13.4e}   {C(5,5):13.4e}')
    store.write(NFEW, f'{C(3,6):13.4e}   {C(5,6):13.4e}    {C(6,6):13.4e}')
    if NFEW > 1:
        store.write(
            NFEW,
            f'{C(3,7):13.4e}   {C(5,7):13.4e}   {C(6,7):13.4e}    {C(7,7):13.4e}')  # noqa E501
        store.write(
            NFEW,
            f'{C(3,8):13.4e}   {C(5,8):13.4e}   {C(6,8):13.4e}   {C(7,8):13.4e}   {C(8,8):13.4e}')  # noqa E501

    if NFEW > 2:
        store.write(
            NFEW,
            f'{C(3,9):13.4e}   {C(5,9):13.4e}   {C(6,9):13.4e}   {C(7,9):13.4e}   {C(8,9):13.4e}   {C(9,9):13.4e}')  # noqa #E501
        store.write(
            NFEW,
            f'{C(3,10):13.4e}   {C(5,10):13.4e}   {C(6,10):13.4e}    {C(7,10):13.4e}     {C(8,10):13.4e}   {C(9,10):13.4e}    {C(10,10):13.4e}')  # noqa E501

    store.write(NFEW, ' -------------------------------------------------')
    store.close(unit=1)
    store.close(unit=2)
    store.close(unit=3)


# **<search for one more & refine amplitudes>****************************
def find_latest_peak(COMS, GRAD, HESS, d_params, INDX, COVAR, o_w1,
                     prog, o_bgd, o_el, store, lptfile, DETLOG,
                     fit_and_chi_func):
    # assume that the previous loops have found good estimates for the other
    # quasi elastic peaks-> only need to focus on getting the latest one
    if o_w1 == 1 and COMS["FIT"].NFEW >= 1:
        COMS["FIT"].FITP.set(5, 0.1)
        COMS["FIT"].FITP.set(
            6,
            COMS["QW1"].QW1(
                COMS["QW1"].ISPEC) /
            COMS["SCL"].WSCL)  # est width of elastic peak
        if COMS["FIT"].NFEW == 1:
            return refine_param_values(GRAD, HESS, 3 + COMS["FIT"].NFEW,
                                       DETLOG, INDX, COVAR, COMS,
                                       fit_and_chi_func, prog, o_bgd,
                                       o_w1, o_el, store, lptfile)

    # 4 default params (BG and elastic line) plus 2 per quasi elastic peak
    J = 4 + 2 * COMS["FIT"].NFEW
    weight = 0.85
    N_iterations = NINT(np.log(5.0 * COMS["SCL"].GSCL /
                                COMS["SCL"].WSCL) /
                        np.log(weight))  # int(5*log(offset/width)/log(0.85))
    CMIN = 1.0E20
    COMS["FIT"].FITP.set(J - 1, 0.1)
    # set the params for the last peak being considered
    COMS["FIT"].FITP.set(J, 1.0)
    tmp_width = COMS["FIT"].FITP(J)

    # "minimize" to find the latest peak (and update the previous peaks)
    for II in get_range(1, N_iterations):
        HESS, COVAR, d_params, DETLOG = refine_param_values(
            GRAD, HESS, 3 + COMS["FIT"].NFEW, DETLOG, INDX, COVAR, COMS,
            fit_and_chi_func, prog, o_bgd, o_w1, o_el, store, lptfile)

        CNORM = fit_and_chi_func(COMS["FIT"].FITP, COMS, o_bgd, o_w1)
        if CNORM < CMIN:
            CMIN = CNORM
            tmp_width = COMS["FIT"].FITP(J)  # only keep the "minimal value"

        # makes a step in the "minimization"
        COMS["FIT"].FITP.set(J, COMS["FIT"].FITP(J) * weight)

    # make sure we use the "best value"
    COMS["FIT"].FITP.set(J, tmp_width)
    return refine_param_values(GRAD, HESS, 3 + COMS["FIT"].NFEW,
                               DETLOG, INDX, COVAR, COMS,
                               fit_and_chi_func, prog, o_bgd,
                               o_w1, o_el, store, lptfile)


@deprecated
def SEARCH(COMS, GRAD, HESS, DPAR, INDX, COVAR, o_w1,
           prog, o_bgd, o_el, store, lptfile, DETLOG,
           Chi_func):
    return find_latest_peak(COMS, GRAD, HESS, DPAR, INDX,
                            COVAR, o_w1, prog, o_bgd, o_el,
                            store, lptfile, DETLOG, Chi_func)
