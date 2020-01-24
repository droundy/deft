#!/usr/bin/python2

import numpy

def e_hist(fbase):
    try:
        trans = numpy.loadtxt(fbase +"-transitions.dat", dtype=numpy.float)
        energy = -trans[:, 0]
        trans = trans[:, 1:]
        hist = numpy.sum(trans, axis=1)
    except:
        # energy histogram file; indexed by [-energy,counts]
        e_hist = numpy.loadtxt(fbase+"-E.dat", ndmin=2, dtype=numpy.float)
        energy = -e_hist[:, 0] # array of energies
        hist = e_hist[:, 1]
    return energy, hist

def e_lndos(f):
    if '.dat' not in f:
        f = f+"-lndos.dat"
    e_lndos = numpy.loadtxt(f, ndmin=2, dtype=numpy.float)
    energy = -e_lndos[:, 0] # array of energies
    lndos = e_lndos[:, 1]
    return energy, lndos

def e_lndos_ps(fbase):
    if '.dat' not in fbase:
        fbase = fbase + "-lndos.dat"
    e_lndos_ps = numpy.loadtxt(fbase, ndmin=2, dtype=numpy.float)
    energy = -e_lndos_ps[:, 0]
    lndos = e_lndos_ps[:, 1]
    ps = e_lndos_ps[:, 2] # pessimistic samples

    return energy, lndos, ps

def e_lndos_ps_lndostm(fbase):
    if '.dat' not in fbase:
        fbase = fbase + "-lndos.dat"
    e_lndos_ps = numpy.loadtxt(fbase, ndmin=2, dtype=numpy.float)
    energy = -e_lndos_ps[:, 0]
    lndos = e_lndos_ps[:, 1]
    ps = e_lndos_ps[:, 2] # pessimistic samples
    lndostm = None
    if len(e_lndos_ps[0,:]) >= 4:
        lndostm = e_lndos_ps[:, 3]

    return energy, lndos, ps, lndostm

def e_lnw(fbase):
    e_lnw = numpy.loadtxt(fbase+"-lnw.dat", ndmin=2, dtype=numpy.float)

    energy = -e_lnw[:, 0] # array of energies
    lnw = e_lnw[:, 1]
    return energy, lnw

def T_u_cv_s_minT(fbase):
    max_T = 20.0
    T_bins = 1e3
    dT = max_T/T_bins
    T_range = numpy.arange(dT, max_T, dT)
    
    # Now compute (or just read in) the lndos and the energies
    try:
        energy, ln_dos = e_lndos(fbase)
        min_T = minT(fbase+'-lndos.dat')
    except:
        min_T = minT(fbase)
        # energy histogram file; indexed by [-energy,counts]
        e_hist = numpy.loadtxt(fbase+"-E.dat", ndmin=2)
        # weight file; indexed by [-energy,ln(weight)]
        lnw_hist = numpy.loadtxt(fbase+"-lnw.dat", ndmin=2)

        energy = -e_hist[:, 0] # array of energies
        lnw = lnw_hist[e_hist[:, 0].astype(int), 1] # look up the lnw for each actual energy
        ln_dos = numpy.log(e_hist[:, 1]) - lnw

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

def minT(f):
    if '.dat' not in f:
        f = f+"-E.dat"
    min_T = 0
    with open(f) as file:
        for line in file:
            if("min_T" in line):
                this_min_T = float(line.split()[-1])
                if this_min_T > min_T:
                    min_T = this_min_T
                break
    return min_T

def moves(fbase):
    moves = 0
    with open(fbase) as file:
        for line in file:
            if("total" in line):
                this_moves = float(line.split()[-1])
                if this_moves > moves:
                    moves = this_moves
                break
    return moves

def minT_from_transitions(fbase):
    return minT(fbase+"-transitions.dat")

def convergedT(f):
    if not '.dat' in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "converged temperature:" in line:
                return float(line.split()[-1])
    return 0

def converged_state(f):
    if not '.dat' in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "converged state:" in line:
                return int(line.split()[-1])
    print(('ERROR FINDING converged_state in', f))

def iterations(f):
    if '.dat' not in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "iterations:" in line:
                return int(line.split()[-1])

def wl_factor(f):
    if '.dat' not in f:
        f = f+"-lndos.dat"
    with open(f) as file:
        for line in file:
            if "# WL Factor: " in line:
                return eval(line.split(': ')[-1])

def dr_g(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "de_g" in line:
                return float(line.split()[-1])

def dimensions(f):
    if '.dat' not in f:
        f = f+"-E.dat"
    with open(f) as file:
        for line in file:
            if "# cell dimensions: " in line:
                return eval(line.split(': ')[-1])

def read_N(f):
    if not '.dat' in f:
        f = f+"-transitions.dat"
    with open(f) as file:
        for line in file:
            if "N" in line:
                return int(line.split(': ')[-1])

def read_ff(fbase):
    with open(fbase+"-E.dat") as file:
        for line in file:
            if "ff" in line:
                return float(line.split(': ')[-1])

def min_important_energy(f):
    if not '.dat' in f:
        f = f+"-transitions.dat"
    with open(f) as file:
        for line in file:
            if("min_important_energy" in line):
                return float(line.split()[-1])

def too_low_high_energy(f):
    if not '.dat' in f:
        f = f+"-transitions.dat"
    too_low = None
    too_high = None
    with open(f) as file:
        for line in file:
            if("too_low_energy" in line or "too_lo_energy" in line):
                too_low = float(line.split()[-1])
                if too_high is not None:
                  return too_low, too_high
            if("too_high_energy" in line or "too_hi_energy" in line):
                too_high = float(line.split()[-1])
                if too_low is not None:
                  return too_low, too_high
    return too_low, too_high

def max_entropy_state(f):
    if not '.dat' in f:
        f = f+"-transitions.dat"
    with open(f) as file:
        for line in file:
            if("max_entropy_state" in line):
                return int(line.split()[-1])

def g_r(fbase, T):
    data = numpy.loadtxt(fbase+"-g.dat", ndmin=2)
    ghist = data[1:, 3:]
    dr = dr_g(fbase)
    r_1d = data[0, 3:]
    E_1d = data[1:, 0]
    hist_1d = data[1:, 1]
    lnw_1d = data[1:, 2]
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
    for i in range(n_r):
        dos_boltz[:, i] /= sum(dos_boltz[:, i])

    dr = r[0, 1] - r[0, 0]
    dV = (4/3)*numpy.pi*((r+dr/2)**3 - (r-dr/2)**3)
    ff = read_ff(fbase)
    n = ff/(4/3*numpy.pi)
    N = read_N(fbase)
    g_of_E = ghist/dV/hist/n/N

    g = numpy.zeros(n_r)

    counts = 0
    for i in range(n_E):
        if hist[i, 0]:
            g += dos_boltz[i,:]*g_of_E[i,:]
            counts += dos_boltz[i, 0]*hist[i, 0]
            # if dos_boltz[i,0] > 0.0001:
            #     print 'hist %d at energy %d with weight %g' % (hist[i,0], E[i,0], dos_boltz[i,0])

    # print 'with N %d we have counts %g' % (N, counts)
    return g, r[0,:]


def density_x(fdensity, T):
    if '.dat' not in fdensity:
        fdensity = fdensity+"-density.dat"
    #fdensity = fbase+"-density.dat"
    data = numpy.loadtxt(fdensity)
    denshist = data[1:, 2:]
    x_1d = data[0, 2:]
    dx = x_1d[1] - x_1d[0]
    E_1d = data[1:, 0]
    ln_dos_1d = data[1:, 1]
    hist_1d = numpy.zeros_like(ln_dos_1d)
    n_E = len(E_1d)
    n_x = len(x_1d)
    N = read_N(fdensity)
    for i in range(n_E):
        hist_1d[i] = sum(denshist[i,:])/N
    x, E = numpy.meshgrid(x_1d, E_1d)
    x, ln_dos = numpy.meshgrid(x_1d, ln_dos_1d)
    x, hist = numpy.meshgrid(x_1d, hist_1d)

    ln_dos_boltz_1d = ln_dos_1d - E_1d/T
    dos_boltz_1d = numpy.exp(ln_dos_boltz_1d - ln_dos_boltz_1d.max())
    dos_boltz_1d[hist_1d == 0] = 0

    # now let's normalize the density of states
    dos_boltz_1d /= sum(dos_boltz_1d)

    # note:  what w
    lenx, leny, lenz = dimensions(fdensity)
    density_of_E = denshist/(dx*leny*lenz*hist)*(4*numpy.pi/3)

    density = numpy.zeros(n_x)

    counts = 0
    for i in range(n_E):
        if hist_1d[i]:
            density += dos_boltz_1d[i]*density_of_E[i,:]
            counts += dos_boltz_1d[i]*hist_1d[i]
            # if dos_boltz_1d[i] > 0.0001:
            #     print 'hist %d at energy %d with weight %g' % (hist_1d[i], E_1d[i], dos_boltz_1d[i])

    # print 'we have counts %g with element dimensions %gx%gx%g' % (counts, dx, leny, lenz)
    # print 'sum', sum(hist_1d), sum(sum(denshist))
    return density, x_1d

def e_de_transitions(basename):
    trans = numpy.loadtxt(basename+"-transitions.dat", dtype=numpy.float)
    with open(basename+"-transitions.dat") as f:
        for line in f:
            if '# energy\t' in line:
                de = numpy.array([ -float(val) for val in line.split()[2:] ])
                break
    N = read_N(basename)
    e = -trans[:, 0]/N
    de /= -N
    trans = trans[:, 1:]
    for i in range(len(e)):
        trans[i,:] /= sum(trans[i,:])
    e, de = numpy.meshgrid(e, de)
    return e, de, trans

def total_init_iterations(basename):
    trans = numpy.loadtxt(basename+"-transitions.dat", dtype=numpy.float)
    trans = trans[:, 1:]
    return numpy.sum(trans)/read_N(basename)

def e_and_total_init_histogram(basename):
    trans = numpy.loadtxt(basename+"-transitions.dat", dtype=numpy.float)
    N = read_N(basename)
    e = -trans[:, 0]
    trans = trans[:, 1:]
    return e, numpy.sum(trans, axis=1)

def e_diffusion_estimate(basename):
    e, de, trans = e_de_transitions(basename)
    e = e[0,:]
    de = de[:, 0]
    diffusion = numpy.zeros_like(e)
    for i in range(len(e)):
        meane = sum(de*trans[i,:])/len(de)
        meane2 = sum(de**2*trans[i,:])/len(de)
        norm = sum(trans[i,:])
        diffusion[i] = numpy.sqrt(meane2 - meane**2)
    return e, diffusion
