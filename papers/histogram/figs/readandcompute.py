#!/usr/bin/python2
import numpy

def u_cv_s(ww, ff, N, method):
    max_T = 1.4
    T_bins = 1e3
    dT = max_T/T_bins
    T_range = numpy.arange(dT,max_T,dT)
    min_T = 0
    try:
        with open("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, method)) as file:
            for line in file:
                if("min_T" in line):
                    this_min_T = float(line.split()[-1])
                    if this_min_T > min_T:
                        min_T = this_min_T
                    break
        # energy histogram file; indexed by [-energy,counts]
        e_hist = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, method), ndmin=2)
        # weight histogram file; indexed by [-energy,ln(weight)]
        lnw_hist = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, method), ndmin=2)
    except:
        # the files must not exist
        return None

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = numpy.log(e_hist[:,1]) - lnw

    Z = numpy.zeros(len(T_range)) # partition function
    U = numpy.zeros(len(T_range)) # internal energy
    CV = numpy.zeros(len(T_range)) # heat capacity
    S = numpy.zeros(len(T_range)) # entropy

    Z_inf = sum(numpy.exp(ln_dos - ln_dos.max()))
    S_inf = sum(-numpy.exp(ln_dos - ln_dos.max())*(-ln_dos.max() - numpy.log(Z_inf))) / Z_inf

    for i in range(len(T_range)):
        ln_dos_boltz = ln_dos - energy/T_range[i]
        dos_boltz = numpy.exp(ln_dos_boltz - ln_dos_boltz.max())
        Z[i] = sum(dos_boltz)
        U[i] = sum(energy*dos_boltz)/Z[i]
        # S = \sum_i^{microstates} P_i \log P_i
        # S = \sum_E D(E) e^{-\beta E} \log\left(\frac{e^{-\beta E}}{\sum_{E'} D(E') e^{-\beta E'}}\right)
        S[i] = sum(-dos_boltz*(-energy/T_range[i] - ln_dos_boltz.max() \
                                       - numpy.log(Z[i])))/Z[i]
        # Actually compute S(T) - S(T=\infty) to deal with the fact
        # that we don't know the actual number of eigenstates:
        S[i] -= S_inf
        CV[i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2
    return U, CV, S
