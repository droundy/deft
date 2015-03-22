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
#arg methods = [["nw","wang_landau","simple_flat","optimized_ensemble","tmmc","oetmmc"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-ps.dat" % (ww, ff, N, method) for method in methods]

plt.title('Iterations per sample for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))
for method in methods:
    data_file = "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-ps.dat" % (ww, ff, N, method)
    with open(data_file,'r') as file_handle:
        for line in file_handle:
            entries = line.split()
            if 'iterations:' in entries:
                iterations = int(entries[entries.index('iterations:')+1].replace(',',''))
                continue
    data = numpy.loadtxt(data_file)
    energy = -data[:,0]/N
    round_trips = data[:,1]
    energy = energy[round_trips != 0]
    round_trips = round_trips[round_trips != 0]
    if sum(round_trips) > 0:
        plt.semilogy(energy, iterations/round_trips, styles.dots(method),
                     label=styles.title(method))

plt.xlabel('$U/N\epsilon$')
plt.ylabel('Iterations per sample')
plt.legend(loc='best').get_frame().set_alpha(0.25)
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-sample-rate.pdf" % (ww*100, ff*100, N))

