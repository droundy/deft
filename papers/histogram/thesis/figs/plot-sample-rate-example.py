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

if len(sys.argv) != 6:
    print 'useage: %s ww ff N seed method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = int(sys.argv[3])
#arg N = [25]

seed = int(sys.argv[4])
#arg seed = [0]

methods = eval(sys.argv[5])
#arg methods = [["nw","cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-ps.dat" % (seed, ww, ff, N, method) for method in methods]

fig = plt.figure('ps',figsize=(6.6,4))
ax = fig.add_axes([0.1, 0.1, 0.55, 0.85])

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
    energy = -data[:,0]
    samples = data[:,1]
    energy = energy[samples != 0]
    samples = samples[samples != 0]
    if sum(samples) > 0:
        plt.semilogy(energy, iterations/samples, styles.dots(method),
                       label=styles.title(method))

plt.xlabel('$E/\epsilon$')
plt.ylabel('$s_p(E)$')

ax.legend(loc='upper right',bbox_to_anchor=(1.65,0.95))
plt.savefig("figs/sample-rate-example.pdf")
