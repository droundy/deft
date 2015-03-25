#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N methods show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

# note: speficic HC should be independent of N, but we have to choose one
N = float(sys.argv[3])
#arg N = range(5,21)

methods = eval(sys.argv[4])
#arg methods = [["tmmc-golden","wang_landau","simple_flat","tmmc","oetmmc","wang_landau_oe","simple_flat_oe","tmmc_oe","oetmmc_oe"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

reference = "tmmc-golden"

max_T = 2
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)
min_T = 0 # we will adjust this

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity
S = {} # entropy

# we want to keep our methods distinct from our reference
if reference in methods:
    methods.remove(reference)

for method in set(methods+[reference]):

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

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = numpy.log(e_hist[:,1]) - lnw

    Z = numpy.zeros(len(T_range)) # partition function
    U[method] = numpy.zeros(len(T_range)) # internal energy
    CV[method] = numpy.zeros(len(T_range)) # heat capacity
    S[method] = numpy.zeros(len(T_range)) # entropy

    Z_inf = sum(numpy.exp(ln_dos - ln_dos.max()))
    S_inf = sum(-numpy.exp(ln_dos - ln_dos.max())*(-ln_dos.max() - numpy.log(Z_inf))) / Z_inf

    for i in range(len(T_range)):
        ln_dos_boltz = ln_dos - energy/T_range[i]
        dos_boltz = numpy.exp(ln_dos_boltz - ln_dos_boltz.max())
        Z[i] = sum(dos_boltz)
        U[method][i] = sum(energy*dos_boltz)/Z[i]
        # S = \sum_i^{microstates} P_i \log P_i
        # S = \sum_E D(E) e^{-\beta E} \log\left(\frac{e^{-\beta E}}{\sum_{E'} D(E') e^{-\beta E'}}\right)
        S[method][i] = sum(-dos_boltz*(-energy/T_range[i] - ln_dos_boltz.max() \
                                       - numpy.log(Z[i])))/Z[i]
        # Actually compute S(T) - S(T=\infty) to deal with the fact
        # that we don't know the actual number of eigenstates:
        S[method][i] -= S_inf
        CV[method][i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2

    plt.figure('u')
    plt.plot(T_range,U[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('hc')
    plt.plot(T_range,CV[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('S')
    plt.plot(T_range,S[method]/N,styles.plot(method),label=styles.title(method))

for method in methods:

    plt.figure('u_err')
    plt.plot(T_range,(U[method]-U[reference])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('hc_err')
    plt.plot(T_range,(CV[method]-CV[reference])/N,
             styles.plot(method),label=styles.title(method))

    plt.figure('S_err')
    plt.plot(T_range,(S[method]-S[reference])/N,
             styles.plot(method),label=styles.title(method))

plt.figure('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))

plt.figure('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))

plt.figure('S')
plt.title('Configurational entropy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-S.pdf" % (ww*100, ff*100, N))

plt.figure('u_err')
plt.title('Error in specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u_err.pdf" % (ww*100, ff*100, N))

plt.figure('hc_err')
plt.title('Error in specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$\\Delta C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc_err.pdf" % (ww*100, ff*100, N))

plt.figure('S_err')
plt.title('Error in configurational entropy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))
plt.xlabel('$kT/\epsilon$')
plt.ylabel(r'$\Delta S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-S_err.pdf" % (ww*100, ff*100, N))

min_T_i = int(min_T/max_T*T_bins)
error_data = open("figs/error-table-ww%02.0f-ff%02.0f-%i.dat" % (ww*100, ff*100, N), "w")
error_data.write("# min_T: %g\n" % min_T)
error_data.write("# %18s %18s %18s %18s\n"
                 % ('method', 'u_error (T)', 'cv_error (T)', 'S_error (T)'))

for method in methods:

    dU = abs(U[method][min_T_i:] - U[reference][min_T_i:])
    dCV = abs(CV[method][min_T_i:] - CV[reference][min_T_i:])
    dS = abs(S[method][min_T_i:] - S[reference][min_T_i:])

    u_error = max(dU)/N
    cv_error = max(dCV)/N
    s_error = max(dS)/N

    u_error_T = T_range[numpy.argmax(dU)]
    cv_error_T = T_range[numpy.argmax(dCV)]
    s_error_T = T_range[numpy.argmax(dS)]

    error_data.write("%20s %10g (%.3f) %10g (%.3f) %10g (%.3f)\n"
                     % (method, u_error, u_error_T, cv_error, cv_error_T, s_error, s_error_T))

error_data.close()

