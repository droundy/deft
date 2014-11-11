#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, glob
import styles

if len(sys.argv) not in [4,5]:
    print 'useage: %s ww ff N versions show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3, 0.4]

versions = eval(sys.argv[3])
#arg versions = [["nw", "wang_landau", "gaussian", "flat", "walkers", "kT2", "kT1"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-rt.dat" % (ww, ff, N, version) for version in versions for N in [20]]

# input: "data/periodic-ww*-ff*-N*-*-rt.dat"

plt.title('Scaling for $\lambda=%g$, $\eta=%g$' % (ww, ff))
for version in versions:
    for filename in glob.glob("data/periodic-ww%04.2f-ff%04.2f-N*-%s-rt.dat" % (ww, ff, version)):
        data = numpy.loadtxt(filename)
        plt.plot([1,2], [3,4], label=filename)

plt.xlabel('$U/N\epsilon$')
plt.ylabel('Round trips')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-scaling.pdf" % (ww*100, ff*100))

if 'show' in sys.argv:
    plt.show()
