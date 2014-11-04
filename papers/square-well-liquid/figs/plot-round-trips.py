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
#arg versions = [["nw", "wang_landau", "gaussian", "flat", "walkers", "kT2", "kT1"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-rt.dat" % (ww, ff, N, version) for version in versions]

plt.title('Round trips for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

for version in versions:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-rt.dat" % (ww, ff, N, version))
    energy = -data[:,0]/N
    round_trips = data[:,1]
    energy = energy[round_trips != 0]
    round_trips = round_trips[round_trips != 0]
    if sum(round_trips) > 0:
        plt.semilogy(energy, round_trips, styles.dots[version], label=styles.title[version])

plt.xlabel('$U/N\epsilon$')
plt.ylabel('Round trips')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-rt.pdf" % (ww*100, ff*100, N))

if 'show' in sys.argv:
    plt.show()
