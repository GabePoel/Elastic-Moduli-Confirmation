# fderivs corresponds to df/dc
# fderivs includes derivatives with respect to the sample dimensions (which my code will not need to do)
# why eight terms, are two of them sample dimensions?

import numpy as np
import RUS_fortran_functions as funcs
from RUS_basics import forwardsolve, basisfuncs, build_C, assign, rotate_C_general


def build_alpha(
        indep_elems,
        gfreqs,
        basis,
        mass,
        density,
        shape,
        dims,
        fitdims):

    # if density is known, use that to calculate mass
    if density > 0:
        [Lx, Ly, Lz] = dims
        if shape == 1:  # RPR
            mass = density * 8 * Lx * Ly * Lz
        if shape == 2:  # diamond prism
            mass = density * 4 * Lx * Ly * Lz
        if shape == 3:  # cylinder
            mass = density * (2 * Lz) * np.pi * Lx**2

    # building the rank-4 tensor C, from the independent elements given
    C = build_C(indep_elems)
    freqs, eigenvects = forwardsolve(
        basis, mass, C, shape, dims)  # MEMORYHOG (more minor)
    print('freqs')
    print(freqs)
    print(freqs.shape)

    # cut down length of freqs so that "assignment" can go fairly quickly
    high_gfreq = gfreqs[-1]
    print('high_gfreq')
    print(high_gfreq)
    print(high_gfreq.shape)

    # index of first calculated frequency greater than 2*highest measured
    # frequency
    double_freq = np.where(freqs > 2 * high_gfreq)[0][0]
    twice_number = 2 * len(gfreqs)  # twice the number of measured frequencies
    # use whichever is greater for the cutoff
    cutoff = max(double_freq, twice_number)
    freqs = freqs[0:cutoff]
    eigenvects = eigenvects[:, 0:cutoff]

    # print len(freqs)
    #pairs = assignment.Hungarian(freqs,gfreqs)

    # "pairs" consists of two columns: the first is just the range of 0 to 1
    # less than the length of gfreqs (indices of measured freqs, plus extra
    # "dummies" where resonances are missed).  the second is the indices of the
    # calculated freqs with which the gfreqs are matched.
    # the value of "pairs[k+1,1]-pairs[k,1] - 1" indicates the number of
    # missing modes between gfreqs[k] and gfreqs[k+1].  (it should be 0 if
    # there are no missing modes.)

    # go through measured frequencies; add dummy frequencies where the
    # assignment algorithm indicates they will improve the fit, and set the
    # weights of these dummy frequencies to 0

    Nfreqs = len(gfreqs)
    # cut down to # of freqs we've measured plus # of dummy freqs inserted
    freqs = freqs[0:Nfreqs]
    eigenvects = eigenvects[:, 0:Nfreqs]

    parameters = indep_elems  # Removed fit dimensions from paramaters
    print(parameters)
    # return 0

    # create a list of dGamma/dp, one for each p
    R = 3 * len(basis)
    Gammaderivs = np.zeros((len(parameters), R, R))  # MEMORYHOG
    Ederivs = np.zeros((len(parameters), R, R))  # MEMORYHOG
    for pp in range(0, len(indep_elems)):  # derivatives w.r.t. elastic moduli
        # make a list of elastic constants that is zero except the one in
        # question, which is 1
        new_elems = np.zeros(len(indep_elems))
        new_elems[pp] = 1
        # build the rank-4 tensor of these 0/1 elastic constants
        Cforderiv = build_C(new_elems)
        # building Gamma from this C is like taking the derivative w.r.t. that
        # elastic constant!
        Gammaderivs[pp], Ederivs[pp] = funcs.matrixderivs(
            basis, mass, Cforderiv, shape, dims, 0, 0)
    # Removed Gammaderivs and Ederivs w.r.t. dimensions

    # build an array of the terms df/dp (each row is for one p, each column is
    # for one f)
    fderivs = np.zeros((len(parameters), len(freqs)))
    print(fderivs)
    for ii in range(0, len(freqs)):
        for pp in range(0, len(parameters)):
            omega = np.multiply(
                (2 * np.pi) / 1000,
                freqs[ii])  # (omega in units of MHz)
            # MEMORYHOG (more minor)
            middle_matrix = Gammaderivs[pp] - omega**2 * Ederivs[pp]
            domega2dp = eigenvects[:, ii].dot(
                middle_matrix).dot(eigenvects[:, ii])
            # (in units of kHz per...parameter change)
            dfdp = domega2dp * (10**6) / (8 * np.pi**2 * freqs[ii])
            #domega2dp = eigenvects[:,ii].dot(Gammaderivs[pp]).dot(eigenvects[:,ii])
            # dfdp = domega2dp/(8*pi**2*freqs[ii])  #(in units of kHz
            # per...parameter change)
            fderivs[pp, ii] = dfdp

    alpha = quick_alpha(fderivs, indep_elems, freqs)

    print('alpha sum')
    for pp in range(0, len(parameters)):
        print('Parameter ' + str(pp) + ':')
        alpha_sum = np.sum(alpha[..., pp])
        print(alpha_sum)

    return alpha


def quick_alpha(fderivs, indep_elems, freqs):
    alpha = np.zeros((len(indep_elems), len(freqs)))
    for ii in range(0, len(freqs)):
        for pp in range(0, len(indep_elems)):
            fderiv_element = fderivs[pp, ii]
            reference_element = indep_elems[pp] / freqs[ii]
            alpha[pp, ii] = fderiv_element * reference_element
    return alpha


def consistency_check():
    # WiP
    pass
