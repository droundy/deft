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
#arg N = [200]

versions = eval(sys.argv[4])
#arg versions = [["-nw", "-flat", "-gaussian", "-kT2", "-kT1", "-kT0.1"]]

max_hc = 3 # upper axis limit when plotting heat capacity

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i%s-%s.dat" % (ww, ff, N, version, data) for version in versions for data in ["E","lnw"]]

# FIXME: make these inputs?
kTmin = 0
kTmax = 10
dkT = float(kTmax-kTmin) / 1000
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

# specific internal energy as a function of kT
def u(kT_array,e_hist,lnw_hist):
    energy = -e_hist[:,0]
    dos = e_hist[:,1]*numpy.exp(-lnw_hist[:,1])
    u_out = numpy.zeros(len(kT_array))
    for i in range(len(u_out)):
        u_out[i] = sum(energy*dos*numpy.exp(-(energy-min(energy))/kT_array[i])) \
          / sum(dos*numpy.exp(-(energy-min(energy))/kT_array[i]))
    return u_out/N

# specific heat capacity as a function of kT
def cv(e_hist,lnw_hist):
    return (u(kT_range+dkT/2,e_hist,lnw_hist)
            - u(kT_range-dkT/2,e_hist,lnw_hist)) / dkT

fig_u, ax_u = sci_fig('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

fig_hc, ax_hc = sci_fig('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

for version in versions:
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i%s-E.dat" % (ww, ff, N, version), ndmin=2)
    lnw_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i%s-lnw.dat" % (ww, ff, N, version), ndmin=2)

    plt.figure('u')
    plt.plot(kT_range,u(kT_range,e_hist,lnw_hist),
             styles.plot[version],label=styles.title[version])

    plt.figure('hc')
    plt.plot(kT_range,cv(e_hist,lnw_hist),styles.plot[version],label=styles.title[version])

plt.figure('u')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
if 'show' in sys.argv:
    plt.show()
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))
plt.close()

plt.figure('hc')
plt.ylim(0,max_hc)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
if 'show' in sys.argv:
    plt.show()
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))
plt.close()
