#!/usr/bin/python2
from __future__ import division
import numpy

def T_u_cv_s_minT(fbase):
    max_T = 100.0
    T_bins = 1e3
    dT = max_T/T_bins
    T_range = numpy.arange(dT,max_T,dT)
    min_T = minT(fbase)
    # energy histogram file; indexed by [-energy,counts]
    e_hist = numpy.loadtxt(fbase+"-E.dat", ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = numpy.loadtxt(fbase+"-lnw.dat", ndmin=2)

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
    return T_range, U, CV, S, min_T

def minT(fbase):
    min_T = 0
    with open(fbase+"-E.dat") as file:
        for line in file:
            if("min_T" in line):
                this_min_T = float(line.split()[-1])
                if this_min_T > min_T:
                    min_T = this_min_T
                break
    return min_T

def g_r(fbase, ff, T, N):
    E = numpy.loadtxt(fbase+"-gE.dat", ndmin=2)
    lnw = numpy.loadtxt(fbase+"-glnw.dat", ndmin=2)
    r = numpy.loadtxt(fbase+"-gr.dat", ndmin=2)
    hist = numpy.loadtxt(fbase+"-gEhist.dat", ndmin=2)
    ghist = numpy.loadtxt(fbase+"-g.dat", ndmin=2)

    ln_dos = numpy.log(hist) - lnw

    n_E = len(E[:,0])
    n_r = len(E[0,:])

    ln_dos_boltz = ln_dos - E/T
    dos_boltz = numpy.exp(ln_dos_boltz - ln_dos_boltz.max())

    # now let's normalize the density of states
    for i in xrange(n_r):
        dos_boltz[:,i] /= sum(dos_boltz[:,i])

    dr = r[0,1] - r[0,0]
    dV = (4/3)*numpy.pi*((r+dr/2)**3 - (r-dr/2)**3)
    n = ff/(4/3*numpy.pi)
    g_of_E = ghist/dV/hist/n

    g = numpy.zeros(n_r)

    for i in xrange(n_E):
        g += dos_boltz[i,:]*g_of_E[i,:]

    return g, r[0,:]
