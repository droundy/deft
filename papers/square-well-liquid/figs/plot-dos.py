#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import styles

if len(sys.argv) != 5:
    print 'useage: %s ww ff N versions' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [200, 1000]

versions = eval(sys.argv[4])
#arg versions = [["-nw", "-flat", "-gaussian", "-kT2", "-kT1", "-kT0.1"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i%s-%s.dat" % (ww, ff, N, version, data) for version in versions for data in ["E","lnw"]]

plt.title('Density of states for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
for version in versions:
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i%s-E.dat" % (ww, ff, N, version))
    lnw_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i%s-lnw.dat" % (ww, ff, N, version))
    energy = -e_hist[:,0]/N
    lnw_mean = lnw_hist[:,1].mean()
    dos = e_hist[:,1]*numpy.exp(-(lnw_hist[:,1] - lnw_mean))
    dos /= sum(dos)
    plt.semilogy(energy,dos,styles.dots[version],label=styles.title[version])

plt.xlabel('$E/N\epsilon$')
plt.ylabel('$DOS$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-dos.pdf" % (ww*100, ff*100, N))
