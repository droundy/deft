#!/usr/bin/python3

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
from pylab import *
from scipy.special import erf
import matplotlib.lines as mlines
import os

import styles

figscale=1
plt.figure(figsize=(9*figscale, 14.5*figscale))

def plot_radial(reduced_density, temps):
    have_labelled_bh = False
    sigma_over_R=2**(5/6)
    have_labelled_dft = False
    for temp in temps:
        fname = 'figs/new-data/radial-wca-%06.4f-%04.2f.dat' % (temp, reduced_density/100.0)
        data = loadtxt(fname)
        r = data[:, 0]
        nreduced_density = data[:, 1]
        g = nreduced_density/(reduced_density/100.0)
        if have_labelled_dft:
            plot(r, g, styles.new_dft_linestyle(), color = styles.color_from_kT(temp))
        else:
            plot(r, g, styles.new_dft_linestyle(), color = styles.color_from_kT(temp), label = 'SFMT $T^* = %g$' % temp)
            have_labelled_dft = True

        fname = 'figs/new-data/radial-bh-wca-%06.4f-%04.2f.dat' % (temp, reduced_density/100.0)
        data = loadtxt(fname)
        r = data[:, 0]
        nreduced_density = data[:, 1]
        g = nreduced_density/(reduced_density/100.0)
        if have_labelled_bh:
            plot(r, g, styles.bh_dft_linestyle(), color = styles.color[temp])
        else:
            plot(r, g, styles.bh_dft_linestyle(), color = styles.color[temp], label = 'BH $T^* = %g$' % temp)
            have_labelled_bh = True
        #plot(r, g, styles.bh_dft(temp))
        
        fname = 'figs/mcfcc-%04.4f-%.4f.dat.gradial' % (reduced_density/100.0, temp)
        g = loadtxt(fname)
        plot(g[:, 0]/sigma_over_R, g[:, 1], styles.mcwca_linestyle(), color = styles.color_from_kT(temp), label = 'WCA MC $T^*$ = %g' % temp)
            
    title('Radial distribution function at $n^* = %g$' % (reduced_density/100), fontsize=20)
    xlabel(r'$r/\sigma$', fontsize=18)
    plt.xticks(fontsize=17)
    ylabel('$g(r)$', fontsize=18)
    plt.yticks(fontsize=17)
    legend(fontsize=18)
    xlim(0.5, 3)
    ylim(0)

# input: ['figs/mcfcc-%04.4f-%.4f.dat.gradial' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# input: ['figs/new-data/radial-wca-%06.4f-%04.2f.dat' % (temp, 0.6) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# input: ['figs/new-data/radial-bh-wca-%06.4f-%04.2f.dat' % (temp, 0.6) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# input: ['figs/mcfcc-%04.4f-%.4f.dat.gradial' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/radial-wca-%06.4f-%04.2f.dat' % (temp, 1.0) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/radial-bh-wca-%06.4f-%04.2f.dat' % (temp, 1.0) for temp in [10.0, 5.0, 2.5]]

subplot(2, 1, 1)
plot_radial(60, [10, 5.0, 2.5, 1.0])

subplot(2, 1, 2)
plot_radial(100, [10.0, 5.0, 2.5])

plt.tight_layout();
savefig('figs/radial-distribution.pdf', bbox_inches=0)
print('figs/radial-distribution.pdf')

