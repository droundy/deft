#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N versions show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

# note: speficic HC should be independent of N, but we have to choose one
N = float(sys.argv[3])
#arg N = [20, 100, 200, 1000]

versions = eval(sys.argv[4])
#arg versions = [["nw","wang_landau","vanilla_wang_landau","robustly_optimistic","gaussian","bubble_suppression","walker_optimization","kT2","kT1"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, version, data) for version in versions for data in ["E","lnw"]]

# FIXME: make these inputs?
kTmin = 0
kTmax = 10
dkT = float(kTmax-kTmin) * 1e-3
kT_range = numpy.arange(kTmin+dkT,kTmax,dkT)

# make figure with axes labeled using scientific notation
def sci_fig(handle):
    fig = plt.figure(handle)
    ax = fig.add_subplot(1,1,1)
    fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2,3))
    fmt.set_scientific(True)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax

fig_u, ax_u = sci_fig('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

fig_hc, ax_hc = sci_fig('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

for version in versions:
    # energy histogram file; indexed by [-energy,counts]
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version), ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version), ndmin=2)

    energy = -e_hist[:,0] # array of energies

    # log of DoS(E)*exp(E/kT); indexed by [energy,temperature]
    ln_dos_boltz = numpy.zeros((len(e_hist),len(kT_range)))
    for i in range(len(e_hist)):
        ln_dos_boltz[i,:] = numpy.log(e_hist[i,1]) - lnw_hist[i,1] - energy[i]/kT_range

    Z = numpy.zeros(len(kT_range)) # partition function
    U = numpy.zeros(len(kT_range)) # internal energy
    CV = numpy.zeros(len(kT_range)) # heat capacity
    for i in range(len(kT_range)):
        dos_boltz = numpy.exp(ln_dos_boltz[:,i] - ln_dos_boltz[:,i].max())
        Z[i] = sum(dos_boltz)
        U[i] = sum(energy*dos_boltz)/Z[i]
        CV[i] = sum((energy/kT_range[i])**2*dos_boltz)/Z[i] - \
          (sum(energy/kT_range[i]*dos_boltz)/Z[i])**2

    plt.figure('u')
    plt.plot(kT_range,U/N,styles.plot[version],label=styles.title[version])

    plt.figure('hc')
    plt.plot(kT_range,CV/N,styles.plot[version],label=styles.title[version])

plt.figure('u')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))

plt.figure('hc')
plt.ylim(0)
plt.xlim(0, 5)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))

if 'show' in sys.argv:
    plt.show()
