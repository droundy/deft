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
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = [10, 20, 100, 200, 1000]

versions = eval(sys.argv[4])
#arg versions = [["nw","wang_landau","robustly_optimistic","gaussian","bubble_suppression","walker_optimization", "kT0.4", "kT0.5", "tmmc"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, version, data) for version in versions for data in ["E","lnw"]]

minlog = 0
for version in versions:
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version))
    lnw = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version))
    energy = -e_hist[:,0]/N
    log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10_dos = numpy.log10(e_hist[:,1]) - log10w
    log10_dos -= log10_dos.max()
    if log10_dos.min() < minlog:
        minlog = log10_dos.min()
    plt.plot(energy, log10_dos, styles.dots(version),label=styles.title(version))

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
plt.title('Density of states for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-dos.pdf" % (ww*100, ff*100, N))

plt.figure() # weight functions
minlog = 0
for version in versions:
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version))
    lnw = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version))
    energy = -e_hist[:,0]/N
    log10w = -lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10w -= log10w.max()
    if log10w.min() < minlog:
        minlog = log10w.min()
    plt.plot(energy, log10w, styles.dots(version),label=styles.title(version))
plt.ylim(minlog, 0)
locs, labels = plt.yticks()
newlabels = [tentothe(n) for n in locs]
plt.yticks(locs, newlabels)
plt.ylim(minlog, 0)

plt.xlabel('$U/N\epsilon$')
plt.ylabel('$w$')
plt.title('Weighting functions for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-weights.pdf" % (ww*100, ff*100, N))


if 'show' in sys.argv:
    plt.show()
