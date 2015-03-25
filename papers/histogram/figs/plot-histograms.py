#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N methods show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

N = int(sys.argv[3])
#arg N = range(5,21)+[100, 200, 1000]

methods = eval(sys.argv[4])
#arg methods = [["nw","kT0.4","kT0.5","wang_landau","simple_flat","tmmc","oetmmc","wang_landau_oe","simple_flat_oe","tmmc_oe","oetmmc_oe"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version) for version in methods]

plt.title('Energy histogram for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

for version in methods:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version))
    energy = -data[:,0]/N
    DS = data[:,1]
    plt.semilogy(energy, DS, styles.dots(version), label=styles.title(version))

plt.xlabel('$U/N\epsilon$')
plt.ylabel('$H$')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-E.pdf" % (ww*100, ff*100, N))

