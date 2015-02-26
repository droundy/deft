#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import styles

if len(sys.argv) not in [6,7]:
    print 'useage: %s ww ff N show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [10, 20, 100, 200, 1000]

method1 = sys.argv[4]
#arg method1 = ['wang_landau']

method2 = sys.argv[5]
#arg method2 = ['tmmc']

bothdata = [0,0]

data1 = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-transitions.dat" % (ww, ff, N, method1))
data2 = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-transitions.dat" % (ww, ff, N, method2))

de_values = numpy.linspace(12,-12.0,25)

for method, data in [(method1, data1), (method2, data2)]:
    plt.figure()
    plt.title('Energy transitions from %s for $\lambda=%g$, $n^*=%g$, and $N=%i$' % (method, ww, ff, N))

    energy = -data[:,0]
    maxde = 0
    minde = 0
    transitions = data[:,1:]
    de, e = numpy.meshgrid(de_values, energy)

    for i in range(len(energy)):
        norm = transitions[i,:].sum()
        if norm > 0:
            transitions[i,:] /= norm
        for j in range(len(transitions[i,:])):
            if transitions[i,j] > 0:
                if de[i,j] < minde:
                    minde = de[i,j]
                if de[i,j] > maxde:
                    maxde = de[i,j]
    c = plt.pcolor(de-0.5, e, transitions, vmin=0, vmax=1)
    plt.colorbar(c)
    plt.xlim(minde,maxde)
    plt.ylim(e.min(), e.max())

    plt.xlabel('$\Delta E/\epsilon$')
    plt.ylabel('$E/\epsilon$')

    if method == method1:
        plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-transitions.pdf" % (ww*100, ff*100, N, method1))
    if method == method2:
        plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-transitions.pdf" % (ww*100, ff*100, N, method2))

plt.figure()

if data1[0,0] < data2[0,0]:
    data1 = data1[data2[0,0] - data1[0,0]:,:]
else:
    data2 = data2[data1[0,0] - data2[0,0]:,:]

energy = -data1[:,0]
de, e = numpy.meshgrid(de_values, energy)

cut2 = data1[:,0].min() - data2[:,0].min()

for i in range(len(data1[:,0])):
    if i < len(data2[:,0]):
        for j in range(len(data1[1,:])):
            data1[i,j] -= data2[i,j]

c = plt.pcolor(de-0.5, e, data1[:,1:], cmap='RdBu', vmin=-0.05, vmax = 0.05)
plt.colorbar(c)
plt.xlim(-8,8)
plt.ylim(e.min(), e.max())
plt.xlabel('$\Delta E/\epsilon$')
plt.ylabel('$E/\epsilon$')

plt.title('TM discrepancy between %s and %s\nfor $\lambda=%g$, $n^*=%g$, and $N=%i$' % (method1, method2, ww, ff, N))

plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-%s-%s-compare-transitions.pdf" % (ww*100, ff*100, N, method1, method2))

if 'show' in sys.argv:
    plt.show()
