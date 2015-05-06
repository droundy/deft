#!/usr/bin/python2
import matplotlib, sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys

sys.path.append('../figs')
import styles

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig_size = (5,4)

if len(sys.argv) != 6:
    print 'useage: %s ww ff N seed method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = int(sys.argv[3])
#arg N = [20]

seed = int(sys.argv[4])
#arg seed = [0]

methods = eval(sys.argv[5])
#arg methods = [["nw","simple_flat","wang_landau","tmmc","oetmmc","cfw"]]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

plt.figure(figsize=fig_size)
for method in methods:

    wildfilename = "../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%%s.dat" \
                   % (seed, ww, ff, N, method)
    e_hist = numpy.loadtxt(wildfilename%'E')
    lnw = numpy.loadtxt(wildfilename%'lnw')
    energy = -e_hist[:,0]/N

    log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10_dos = numpy.log10(e_hist[:,1]) - log10w
    log10_dos -= log10_dos.max()

    minlog = 0
    if log10_dos.min() < minlog:
        minlog = log10_dos.min()

    plt.plot(energy, log10_dos, styles.dots(method),label=styles.title(method))

def tentothe(n):
    if n == 0:
        return '$1$'
    if n == 10:
        return '$10$'
    if int(n) == n:
        return r'$10^{%d}$' % n
    return r'$10^{%g}$' % n

plt.ylim(minlog, 0)

ylocs, ylabels = plt.yticks()
newylabels = [tentothe(n) for n in ylocs]
plt.yticks(ylocs, newylabels)

plt.xlabel('$E/N\epsilon$')
plt.ylabel('$\\tilde D(E)$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/dos-sample.pdf")
