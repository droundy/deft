#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import styles

if len(sys.argv) not in [6,7]:
    print 'useage: %s ww ff N min_T methods show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = float(sys.argv[3])
#arg N = range(5,21)

min_T = eval(sys.argv[4])
#arg min_T = [0.1]

methods = eval(sys.argv[5])
#arg methods = [["wang_landau","simple_flat","tmmc","oetmmc"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-conv_T%g-%s.dat" % (ww, ff, N, method, min_T, data) for method in methods for data in ["E","lnw"]]

max_T = 2
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity
S = {} # entropy

minlog = 0
for method in methods:

    e_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-conv_T%g-E.dat"
                           % (ww, ff, N, method, min_T), ndmin=2)
    lnw_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-conv_T%g-lnw.dat"
                             % (ww, ff, N, method, min_T), ndmin=2)

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = numpy.log(e_hist[:,1]) - lnw

    log10w = lnw_hist[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10_dos = numpy.log10(e_hist[:,1]) - log10w
    log10_dos -= log10_dos.max()
    if log10_dos.min() < minlog:
        minlog = log10_dos.min()
    plt.figure('dos')
    plt.plot(energy, log10_dos, styles.dots(method),label=styles.title(method))

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
        S[method][i] = sum(-dos_boltz*(-energy/T_range[i] - ln_dos_boltz.max() \
                                       - numpy.log(Z[i])))/Z[i]
        S[method][i] -= S_inf
        CV[method][i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2

    plt.figure('u')
    plt.plot(T_range,U[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('hc')
    plt.plot(T_range,CV[method]/N,styles.plot(method),label=styles.title(method))

    plt.figure('S')
    plt.plot(T_range,S[method]/N,styles.plot(method),label=styles.title(method))


plt.figure('dos')
plt.ylim(minlog, 0)
locs, labels = plt.yticks()
def tentothe(n):
    if n == 0:
        return '1'
    if n == 10:
        return '10'
    if int(n) == n:
        return r'$10^{%d}$' % n
    return r'$10^{%g}$' % n
newlabels = [tentothe(n) for n in locs]
plt.yticks(locs, newlabels)
plt.ylim(minlog, 0)

plt.xlabel('$U/N\epsilon$')
plt.ylabel('$DoS$')
plt.title('Density of states for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          ' ($kT_{min}/\epsilon=%g$)' % (ww, ff, N, min_T))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-dos-conv-T%02.0f.pdf"
            % (ww*100, ff*100, N, min_T*100))

plt.figure('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          ' ($kT_{min}/\epsilon=%g$)' % (ww, ff, N, min_T))
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u-conv-T%02.0f.pdf"
            % (ww*100, ff*100, N, min_T*100))

plt.figure('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          ' ($kT_{min}/\epsilon=%g$)' % (ww, ff, N, min_T))
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc-conv-T%02.0f.pdf"
            % (ww*100, ff*100, N, min_T*100))

plt.figure('S')
plt.title('Configurational entropy for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          ' ($kT_{min}/\epsilon=%g$)' % (ww, ff, N, min_T))
plt.xlabel(r'$kT/\epsilon$')
plt.ylabel(r'$S_{\textit{config}}/Nk$')
plt.legend(loc='best')
plt.axvline(min_T,linewidth=1,color='k',linestyle=':')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-S-conv-T%02.0f.pdf"
            % (ww*100, ff*100, N, min_T*100))
