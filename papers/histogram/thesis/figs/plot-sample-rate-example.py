#!/usr/bin/env python2
import matplotlib, sys, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys

sys.path.append('../figs')
import styles

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

if len(sys.argv) != 6:
    print 'useage: %s ww ff N seeds method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = int(sys.argv[3])
#arg N = [25]

seeds = eval(sys.argv[4])
#arg seeds = [range(30)]

methods = eval(sys.argv[5])
#arg methods = [["nw","cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

fig = plt.figure('ps',figsize=(6.6,4))
ax = fig.add_axes([0.1, 0.1, 0.55, 0.85])

energy = {}
net_rate = {}

for method in methods:
    energies = []
    rates = []
    for seed in seeds:
        basedir = "../data/"
        fname = "periodic-ww%04.2f-ff%04.2f-N%i-%s-ps.dat" % (ww, ff, N, method)
        if method == 'nw':
            basename = basedir+"s000/"+fname
        else:
            basename = basedir+"s%03d/"%seed+fname

        with open(basename,'r') as file_handle:
            for line in file_handle:
                entries = line.split()
                if 'iterations:' in entries:
                    iterations = int(entries[entries.index('iterations:')+1].replace(',',''))
                    continue
        data = numpy.loadtxt(basename)
        energies.append(-data[:,0])
        samples = data[:,1]

        rates.append(numpy.zeros(len(samples)))
        rates[-1][samples != 0] = iterations/samples[samples != 0]


    min_energy = 0
    for i in range(len(energies)):
        min_e = min(energies[i])
        if min_e < min_energy:
            min_energy = int(min_e)
    energy_num = -min_energy + 1

    energy[method] = -numpy.array(range(energy_num))
    net_rate[method] = numpy.zeros(energy_num)

    for i in range(len(seeds)):
        for j in range(len(rates[i])):
            net_rate[method][-energies[i][j]] += rates[i][j]
    net_rate[method] /= len(seeds)

    energy[method] = energy[method][net_rate[method] != 0]
    net_rate[method] = net_rate[method][net_rate[method] != 0]

for method in methods:
    plt.semilogy(energy[method], net_rate[method], styles.dots(method),
                 markerfacecolor='none', markeredgecolor=styles.color(method),
                 label=styles.title(method))

plt.xlabel('$E/\epsilon$')
plt.ylabel('$s_p(E)$')
ax.legend(loc='upper right',bbox_to_anchor=(1.65,0.95))
plt.savefig("figs/sample-rate-example.pdf")
