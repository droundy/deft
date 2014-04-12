#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 5:
    print 'useage: %s ww ff N kTs' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [200, 1000]

kTs = eval(sys.argv[4])
#arg kTs = [[0.1, 1, 2]]

plt.title('Density of states with well width %g and filling fraction %g' % (ww, ff))

data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-dos.dat" % (ww, ff, N))
energy = -data[:,0]/N
DS = data[:,1]
DS /= sum(DS)
plt.semilogy(energy,DS,'.', label='$kT = \infty$ and $N=%d$' % N)

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT) for kT in kTs]
for kT in kTs:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT))
    energy = -data[:,0]/N
    DS = data[:,1]
    DS /= sum(DS)
    plt.semilogy(energy,DS,'.', label=r'$kT = %g\epsilon$ and $N=%d$' % (kT, N))

plt.xlabel('$E/N\epsilon$')
plt.ylabel('$DOS$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-dos.pdf" % (ww*100, ff*100, N))
