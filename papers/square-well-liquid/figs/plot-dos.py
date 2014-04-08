#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 4:
    print 'useage: %s ww ff N' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [200, 1000]

plt.title('Histogram with well width %g and filling fraction %g' % (ww, ff))

data = numpy.loadtxt("data/periodic-ww%03.1f-ff%04.2f-N%i-nw-E.dat" % (ww, ff, N))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label='$kT = \infty$ and $N=%d$' % N)

data = numpy.loadtxt("data/periodic-ww%03.1f-ff%04.2f-N%i-kT%g-E.dat" % (ww, ff, N, 2))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label=r'$kT = %g\epsilon$ and $N=%d$' % (2, N))


data = numpy.loadtxt("data/periodic-ww%03.1f-ff%04.2f-N%i-kT%g-E.dat" % (ww, ff, N, 1))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label=r'$kT = %g\epsilon$ and $N=%d$' % (1, N))

data = numpy.loadtxt("data/periodic-ww%03.1f-ff%04.2f-N%i-kT%g-E.dat" % (ww, ff, N, 0.1))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label=r'$kT = %g\epsilon$ and $N=%d$' % (0.1, N))


plt.xlabel('$E/N\epsilon$')
plt.ylabel('$D$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-E.pdf" % (ww*10, ff*100, N))

