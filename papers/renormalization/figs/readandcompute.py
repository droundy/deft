#!/usr/bin/python2
import numpy as np
import math
import string
import glob

def T_u_F_cv_s_minT(fbase):
    max_T = 1.4
    T_bins = 1e3
    dT = max_T/T_bins
    m = .25   # Need an appropriate scale for this - look at size of epsilon and sigma

    T_range = np.arange(dT,max_T,dT)
    min_T = minT(fbase)
    # energy histogram file; indexed by [-energy,counts]
    e_hist = np.loadtxt(fbase+"-E.dat", ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = np.loadtxt(fbase+"-lnw.dat", ndmin=2)

    V = dimensions(fbase)[0]**3 # assuming it is cubic
    N = read_N(fbase)

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = np.log(e_hist[:,1]) - lnw

    Z = np.zeros(len(T_range)) # partition function
    U = np.zeros(len(T_range)) # internal energy
    CV = np.zeros(len(T_range)) # heat capacity
    S = np.zeros(len(T_range)) # entropy
    F = np.zeros(len(T_range)) # Free energy

    Z_inf = sum(np.exp(ln_dos - ln_dos.max()))
    S_inf = sum(-np.exp(ln_dos - ln_dos.max())*(-ln_dos.max() - np.log(Z_inf))) / Z_inf
    f_abs = absolute_f(fbase)

    for i in range(len(T_range)):
        ln_dos_boltz = ln_dos - energy/T_range[i]
        dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max())
        Z[i] = sum(dos_boltz)
        U[i] = sum(energy*dos_boltz)/Z[i]

        F[i] = T_range[i]*np.log(Z[i]) - f_abs
        # S = \sum_i^{microstates} P_i \log P_i
        # S = \sum_E D(E) e^{-\beta E} \log\left(\frac{e^{-\beta E}}{\sum_{E'} D(E') e^{-\beta E'}}\right)
        S[i] = sum(-dos_boltz*(-energy/T_range[i] - ln_dos_boltz.max() \
                                       - np.log(Z[i])))/Z[i] \
                                       + np.log((V/N)*(m*T_range[i]/(2*np.pi))**(.5)) + 5/2.0
        # Actually compute S(T) - S(T=\infty) to deal with the fact
        # that we don't know the actual number of eigenstates:
        S[i] -= S_inf
        CV[i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2
    return T_range, U, F, CV, S, min_T


def u_F_cv_s_minT(fbase, T):
    m = .25   # Need an appropriate scale for this - look at size of epsilon and sigma

    min_T = minT(fbase)
    # energy histogram file; indexed by [-energy,counts]
    e_hist = np.loadtxt(fbase+"-E.dat", ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = np.loadtxt(fbase+"-lnw.dat", ndmin=2)

    V = dimensions(fbase)[0]**3 # assuming it is cubic
    N = read_N(fbase)

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = np.log(e_hist[:,1]) - lnw

    Z_inf = sum(np.exp(ln_dos - ln_dos.max()))
    S_inf = sum(-np.exp(ln_dos - ln_dos.max())*(-ln_dos.max() - np.log(Z_inf))) / Z_inf
    f_abs = absolute_f(fbase)

    ln_dos_boltz = ln_dos - energy/T
    dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max())
    Z = sum(dos_boltz)
    U = sum(energy*dos_boltz)/Z

    F = T*np.log(Z) - f_abs
    # S = \sum_i^{microstates} P_i \log P_i
    # S = \sum_E D(E) e^{-\beta E} \log\left(\frac{e^{-\beta E}}{\sum_{E'} D(E') e^{-\beta E'}}\right)
    S = sum(-dos_boltz*(-energy/T - ln_dos_boltz.max() \
                           - np.log(Z)))/Z \
        + np.log((V/N)*(m*T/(2*np.pi))**(.5)) + 5/2.0
    # Actually compute S(T) - S(T=\infty) to deal with the fact
    # that we don't know the actual number of eigenstates:
    S -= S_inf
    CV = sum((energy/T)**2*dos_boltz)/Z - (sum(energy/T*dos_boltz)/Z)**2
    return U, F, CV, S, min_T

def absolute_f(fbase):
    # find the partition function yielding the absolute free energy  using 'absolute/' data
    fbase = fbase[:-4]+ '/absolute/'
    num_files = len(glob.glob(fbase+'*.dat')) # figure out how many files there are
    successes = 0
    total = 0
    ratios = np.zeros(num_files)
    absolute_f = 0
    
   
    print("Num_files is: %s" % num_files)
    for j in xrange(0,num_files):
        if j==0:
            filename = fbase+'%05d-g' % (j)
        else:
            filename = fbase + '%05d' % (j)
        print 'filename is "%s" and j is %d' % (filename, j)
        with open(filename+".dat") as file:
            for line in file:
                if ("N: " in line):
                    N = int(line.split()[-1])
                if("ff_small: " in line):
                    ff = float(line.split()[-1])
                if("failed small checks: " in line):
                    successes = float(line.split()[-1])
                if("valid small checks: " in line):
                    total = float(line.split()[-1])
                    break
        ratios[j] = successes/total
        absolute_f += -np.log(ratios[j])/N

    print("Ratios array is: %s" % ratios)
    print("Calculated absolute_f is: %g" % absolute_f)
    #print("Compare with: %g" % ((4*ff - 3*ff**2)/(1-ff)**2))
    return absolute_f

	
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

def dimensions(f):
    if not '.dat' in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "# cell dimensions: " in line:
                return eval(line.split(': ')[-1])

def read_N(f):
    if not '.dat' in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "N" in line:
                return int(line.split(': ')[-1])
