#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import styles

if len(sys.argv) not in [4,5]:
    print 'useage: %s ww ff N show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [10, 20, 100, 200, 1000]

plt.title('Energy transitions for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-wang_landau-transitions.dat" % (ww, ff, N))
energy = -data[:,0]/N
emax = 0
emin = len(energy)
maxde = 0
minde = 0
transitions = data[:,1:]
de, e = numpy.meshgrid(numpy.linspace(-12,12.0,25), numpy.linspace(0, len(energy)-1, len(energy)))
for i in range(len(energy)):
    norm = transitions[i,:].sum()
    if norm > 0:
        if i > emax:
            emax = i
        if i < emin:
            emin = i
        transitions[i,:] /= norm
    for j in range(len(transitions[i,:])):
        if transitions[i,j] > 0:
            if de[i,j] < minde:
                minde = de[i,j]
            if de[i,j] > maxde:
                maxde = de[i,j]
c = plt.pcolor(de-0.5, -e, transitions)
plt.colorbar(c)
plt.xlim(minde,maxde)
plt.ylim(-emin, -emax)

plt.xlabel('$\Delta E/\epsilon$')
plt.ylabel('$E/\epsilon$')

plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-transitions.pdf" % (ww*100, ff*100, N))

if 'show' in sys.argv:
    plt.show()
