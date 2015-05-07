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

fig_size = numpy.array([5,4])
panel_size_ratio = 1.75

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

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

plt.figure('dos',figsize=fig_size)

all_figs, ((ax_hist, ax_lnw),(ax_dos, ax_legend)) \
          = plt.subplots(2,2,figsize=panel_size_ratio*fig_size)

for method in methods:

    wildfilename = "../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%%s.dat" \
                   % (seed, ww, ff, N, method)
    e_hist = numpy.loadtxt(wildfilename%'E')
    lnw = numpy.loadtxt(wildfilename%'lnw')
    energy = -e_hist[:,0]

    log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10_dos = numpy.log10(e_hist[:,1]) - log10w
    log10_dos -= log10_dos.max()

    minlog = 0
    if log10_dos.min() < minlog:
        minlog = log10_dos.min()

    plt.figure('dos')
    plt.plot(energy, log10_dos, styles.dots(method),
             markerfacecolor='none', markeredgecolor=styles.color(method),
             label=styles.title(method))
    ax_dos.plot(energy, log10_dos, styles.dots(method),
                markerfacecolor='none', markeredgecolor=styles.color(method),
                label=styles.title(method))

    ax_hist.semilogy(energy, e_hist[:,1], styles.dots(method),
                     markerfacecolor='none', markeredgecolor=styles.color(method),
                     label=styles.title(method))

    ax_lnw.plot(energy, log10w, styles.dots(method),
                     markerfacecolor='none', markeredgecolor=styles.color(method),
                     label=styles.title(method))

def tentothe(n):
    if int(n) == n:
        return r'$10^{%d}$' % n
    return r'$10^{%g}$' % n

### density of states
plt.figure('dos')

plt.ylim(minlog, 0)
ylocs, ylabels = plt.yticks()
newylabels = [tentothe(n) for n in ylocs]
plt.yticks(ylocs, newylabels)

plt.xlabel('$E/\epsilon$')
plt.ylabel('$\\tilde D(E)$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/dos-example.pdf")

### all figures
all_figs.canvas.draw()

for ax in [ ax_lnw, ax_dos, ax_hist ]:
    ax.set_xlabel('$E/\epsilon$')
ax_hist.set_ylabel('$H(E)$')

ax_dos.set_ylim(minlog, 0)
for ax in [ ax_lnw, ax_dos ]:
    ylabels = [ int(item.get_text().encode('ascii')[1:-1])
                for item in ax.get_yticklabels() ]
    newylabels = [tentothe(n) for n in ylabels]
    ax.set_yticklabels(newylabels)

ax_lnw.set_ylabel('$w(E)$')
ax_dos.set_ylabel('$\\tilde D(E)$')

ax_legend.axis('off')
handles, labels = ax_dos.get_legend_handles_labels()
ax_legend.legend(handles, labels, loc='center')


all_figs.tight_layout(pad=0.2)
all_figs.savefig('figs/array-example.pdf')
