#!/usr/bin/env python2
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

if len(sys.argv) != 7:
    print 'useage: %s ww ff N min_e seed method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]

ff = float(sys.argv[2])
#arg ff = [0.3]

N = int(sys.argv[3])
#arg N = [25]

# fixme: read from file, instead of taking as input
min_e = int(sys.argv[4])
#arg min_e = [-123]

seed = int(sys.argv[5])
#arg seed = [0]

methods = eval(sys.argv[6])
#arg methods = [["nw","cfw","simple_flat","vanilla_wang_landau","wang_landau","tmmc","oetmmc"]]

# input: ["../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (seed, ww, ff, N, method, data) for method in methods for data in ["E","lnw"]]

plt.figure('dos-thesis',figsize=fig_size)
plt.figure('dos-poster',figsize=fig_size)

ax_placement = [0.11, 0.1, 0.55, 0.85]
fig_size_alt = (6.6,4)

fig_dos = plt.figure('dos-example',figsize=fig_size_alt)
ax_dos = fig_dos.add_axes(ax_placement)
fig_hist = plt.figure('hist-example',figsize=fig_size_alt)
ax_hist = fig_hist.add_axes(ax_placement)
fig_lnw = plt.figure('lnw-example',figsize=fig_size_alt)
ax_lnw = fig_lnw.add_axes(ax_placement)

all_figs, ((ax_hist_all, ax_lnw_all),(ax_dos_all, ax_legend_all)) \
          = plt.subplots(2,2,figsize=panel_size_ratio*fig_size)

cap = 60

def get_arrays(wildfilename):

    e_hist = numpy.loadtxt(wildfilename%'E')
    lnw = numpy.loadtxt(wildfilename%'lnw')
    energy = -e_hist[:,0]

    log10w = lnw[e_hist[:,0].astype(int),1]*numpy.log10(numpy.exp(1))
    log10_dos = numpy.log10(e_hist[:,1]) - log10w
    log10_dos -= log10_dos.max()

    return energy, e_hist, log10w, log10_dos


for method in methods:

    wildfilename = "../data/s%03d/periodic-ww%04.2f-ff%04.2f-N%i-%s-%%s.dat" \
                   % (seed, ww, ff, N, method)

    energy, e_hist, log10w, log10_dos = get_arrays(wildfilename)

    if method != 'wang_landau':
        plt.figure('dos-poster')
        plt.plot(energy, log10_dos, styles.dots(method),
                 markerfacecolor='none', markeredgecolor=styles.color(method),
                 label=styles.title(method))

    for ax in [ ax_dos, ax_dos_all ]:
        ax.plot(energy[log10_dos > -cap], log10_dos[log10_dos > -cap],
                styles.dots(method), markerfacecolor='none',
                markeredgecolor=styles.color(method), label=styles.title(method))

    for ax in [ ax_hist, ax_hist_all ]:
        ax.semilogy(energy, e_hist[:,1], styles.dots(method),
                    markerfacecolor='none', markeredgecolor=styles.color(method),
                    label=styles.title(method))

    for ax in [ ax_lnw, ax_lnw_all ]:
        ax.plot(energy[log10w < cap], log10w[log10w < cap], styles.dots(method),
                markerfacecolor='none', markeredgecolor=styles.color(method),
                label=styles.title(method))

# input: ["../data/periodic-ww1.30-ff0.10-N100-tmmc-golden-%s.dat" % data for data in ["E","lnw"]]

wildfilename = "../data/periodic-ww1.30-ff0.10-N100-tmmc-golden-%s.dat"
energy, e_hist, log10w, log10_dos = get_arrays(wildfilename)

plt.figure('dos-thesis')
plt.plot(energy, log10_dos,'k.')

def tentothe(n):
    if int(n) == n:
        return r'$10^{%d}$' % n
    return r'$10^{%g}$' % n

### density of states for thesis and poster
for dos_plot in ['dos-thesis','dos-poster']:
    plt.figure(dos_plot)

    ylocs, ylabels = plt.yticks()
    newylabels = [tentothe(n) for n in ylocs]
    plt.yticks(ylocs, newylabels)

    plt.xlabel('$E/\epsilon$')
    plt.ylabel('$\\tilde D(E)$')
    plt.tight_layout(pad=0.2)

plt.figure('dos-thesis')
plt.savefig("figs/dos-thesis-example.pdf")

plt.figure('dos-poster')
plt.legend(loc='best')
plt.savefig("figs/dos-poster-example.pdf")

### all figures
all_figs.canvas.draw()

for ax in [ ax_dos, ax_dos_all, ax_hist, ax_hist_all, ax_lnw, ax_lnw_all ]:
    ax.axvline(min_e,linewidth=1,color='k',linestyle=':')
    ax.set_xlabel('$E/\epsilon$')

for ax in [ ax_dos, ax_lnw, ax_dos_all, ax_lnw_all ]:
    newylabels = [tentothe(n) for n in ax.get_yticks()]
    ax.set_yticklabels(newylabels)

for ax in [ ax_dos, ax_dos_all ]:
    ax.set_ylabel('$\\tilde D(E)$')

for ax in [ ax_hist, ax_hist_all ]:
    ax.set_ylabel('$H(E)$')

for ax in [ ax_lnw, ax_lnw_all ]:
    ax.set_ylabel('$w(E)$')

ax_legend_all.axis('off')
handles, labels = ax_dos_all.get_legend_handles_labels()
ax_legend_all.legend(handles, labels, loc='center')

all_figs.tight_layout(pad=0.2)
all_figs.savefig('figs/array-example.pdf')

legend_loc = (1.63,0.8)

plt.figure('dos-example')
ax_dos.legend(loc='upper right',bbox_to_anchor=legend_loc)
plt.savefig('figs/dos-example.pdf')

plt.figure('hist-example')
ax_hist.legend(loc='upper right',bbox_to_anchor=legend_loc)
plt.savefig('figs/hist-example.pdf')

plt.figure('lnw-example')
ax_lnw.legend(loc='upper right',bbox_to_anchor=legend_loc)
plt.savefig('figs/lnw-example.pdf')

