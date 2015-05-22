#!/usr/bin/python2
from __future__ import division
import numpy

def T_u_cv_s_minT(fbase):
    max_T = 20.0
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

def dr_g(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "de_g" in line:
                return float(line.split()[-1])

def dimensions(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "# cell dimensions: " in line:
                return eval(line.split(': ')[-1])

def read_N(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "N" in line:
                return int(line.split(': ')[-1])

def read_ff(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "ff" in line:
                return float(line.split(': ')[-1])

def g_r(fbase, T):
    data = numpy.loadtxt(fbase+"-g.dat", ndmin=2)
    ghist = data[1:,3:]
    dr = dr_g(fbase)
    r_1d = data[0,3:]
    E_1d = data[1:,0]
    hist_1d = data[1:,1]
    lnw_1d = data[1:,2]
    n_E = len(E_1d)
    n_r = len(r_1d)
    r, E = numpy.meshgrid(r_1d, E_1d)
    r, lnw = numpy.meshgrid(r_1d, lnw_1d)
    r, hist = numpy.meshgrid(r_1d, hist_1d)

    ln_dos = numpy.log(hist) - lnw


    ln_dos_boltz = ln_dos - E/T
    dos_boltz = numpy.exp(ln_dos_boltz - ln_dos_boltz.max())
    dos_boltz[hist == 0] = 0

    # now let's normalize the density of states
    for i in xrange(n_r):
        dos_boltz[:,i] /= sum(dos_boltz[:,i])

    dr = r[0,1] - r[0,0]
    dV = (4/3)*numpy.pi*((r+dr/2)**3 - (r-dr/2)**3)
    ff = read_ff(fbase)
    n = ff/(4/3*numpy.pi)
    N = read_N(fbase)
    g_of_E = ghist/dV/hist/n/N

    g = numpy.zeros(n_r)

    counts = 0
    for i in xrange(n_E):
        if hist[i,0]:
            g += dos_boltz[i,:]*g_of_E[i,:]
            counts += dos_boltz[i,0]*hist[i,0]
            if dos_boltz[i,0] > 0.0001:
                print 'hist %d at energy %d with weight %g' % (hist[i,0], E[i,0], dos_boltz[i,0])

    print 'with N %d we have counts %g' % (N, counts)
    return g, r[0,:]


def density_x(fbase, T):
    data = numpy.loadtxt(fbase+"-density.dat")
    denshist = data[1:,2:]
    x_1d = data[0,2:]
    dx = x_1d[1] - x_1d[0]
    E_1d = data[1:,0]
    lnw_1d = data[1:,1]
    hist_1d = numpy.zeros_like(lnw_1d)
    n_E = len(E_1d)
    n_x = len(x_1d)
    N = read_N(fbase)
    for i in xrange(n_E):
        hist_1d[i] = sum(denshist[i,:])/N
    x, E = numpy.meshgrid(x_1d, E_1d)
    x, lnw = numpy.meshgrid(x_1d, lnw_1d)
    x, hist = numpy.meshgrid(x_1d, hist_1d)

    ln_dos_1d = numpy.log(hist_1d) - lnw_1d

    ln_dos_boltz_1d = ln_dos_1d - E_1d/T
    dos_boltz_1d = numpy.exp(ln_dos_boltz_1d - ln_dos_boltz_1d.max())
    dos_boltz_1d[hist_1d == 0] = 0

    # now let's normalize the density of states
    dos_boltz_1d /= sum(dos_boltz_1d)

    # note:  what w
    lenx, leny, lenz = dimensions(fbase)
    density_of_E = denshist/(dx*leny*lenz*hist)*(4*numpy.pi/3)

    density = numpy.zeros(n_x)

    counts = 0
    for i in xrange(n_E):
        if hist_1d[i]:
            density += dos_boltz_1d[i]*density_of_E[i,:]
            counts += dos_boltz_1d[i]*hist_1d[i]
            if dos_boltz_1d[i] > 0.0001:
                print 'hist %d at energy %d with weight %g' % (hist_1d[i], E_1d[i], dos_boltz_1d[i])

    print 'we have counts %g with element dimensions %gx%gx%g' % (counts, dx, leny, lenz)
    print 'sum', sum(hist_1d), sum(sum(denshist))
    return density, x_1d
