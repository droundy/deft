#!/usr/bin/python2
import matplotlib, sys
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys

sys.path.append('../figs')
import styles

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

fig_size = numpy.array([6,5])

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
#arg methods = [["nw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc","cfw"]]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["ps","os"]]

plt.figure('ps',figsize=fig_size)
ax_ps = plt.subplot(111)

for method in methods:
    data_file = "../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-ps.dat" \
                % (seed, ww, ff, N, method)
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

plt.xlabel('$E/N\epsilon$')
plt.ylabel('Iterations per pessimistic sample')
plt.legend(loc='upper right')
plt.tight_layout(pad=0.2)
ax_ps.legend(bbox_to_anchor=(1.16,1))
plt.savefig("figs/sample-rate-sample.pdf")

plt.show()
