#!/usr/bin/python2
import numpy
import math
import string

def T_u_cv_s_minT(fbase):
    max_T = 1.4
    T_bins = 1e3
    dT = max_T/T_bins
    m = .25   # Need an appropriate scale for this - look at size of epsilon and sigma
    V = 1000.0
    N = 10.0

    T_range = numpy.arange(dT,max_T,dT)
    min_T = minT(fbase)
    # energy histogram file; indexed by [-energy,counts]
    e_hist = numpy.loadtxt(fbase+"-E.dat", ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = numpy.loadtxt(fbase+"-lnw.dat", ndmin=2)
      

    with open(fbase+"-E.dat") as file:
        for line in file:
            line = line.translate(cleanse)
            if("dimensions" in line):
                V = line.split()[1]
                print V
            if(" N: " in line):
                N = float(line.split()[-1])
                print N
            break
        

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
                                       - numpy.log(Z[i])))/Z[i] \
                                       + numpy.log((V/N)*(m*T_range[i]/(2*numpy.pi))**(.5)) + 5/2.0
        # Actually compute S(T) - S(T=\infty) to deal with the fact
        # that we don't know the actual number of eigenstates:
        S[i] -= S_inf
        CV[i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2
   # return T_range, U, CV, S, min_T

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
