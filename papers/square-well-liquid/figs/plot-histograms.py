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

plt.title('Histogram for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i%s-E.dat" % (ww, ff, N, version) for version in versions]
for version in versions:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i%s-E.dat" % (ww, ff, N, version))
    energy = -data[:,0]/N
    DS = data[:,1]
    DS /= sum(DS)
    plt.semilogy(energy,DS, styles.dots[version], label=styles.title[version])

plt.xlabel('$E/N\epsilon$')
plt.ylabel('$D$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-E.pdf" % (ww*100, ff*100, N))

