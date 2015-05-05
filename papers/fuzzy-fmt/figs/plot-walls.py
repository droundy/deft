#!/usr/bin/python

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

def plot_walls(reduced_density, temps):
    figure()
    sigma_over_R=2**(5/6)
    have_labelled_dft = False
    for temp in temps:
        fname = 'figs/new-data/wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:,0]
        nreduced_density = data[:,1]
        if have_labelled_dft:
            plot(z, nreduced_density, styles.new_dft_code(temp))
        else:
            plot(z, nreduced_density, styles.new_dft_code(temp), label = 'DFT $T^* = %g$' % temp)
            have_labelled_dft = True
        
        fname = 'figs/mcfcc-walls-%04.4f-%.4f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        zmin = data[:,0].min()
        plot((data[:,0]-zmin)/sigma_over_R, data[:,1]*sigma_over_R**3,
             styles.mcwca(temp), label = 'WCA MC $T^*$ = %g' % temp)

    title('Hard walls with bulk $n^* = %g$' % (reduced_density/100))
    xlabel(r'$z/\sigma$')
    ylabel('$n^*(r)$')
    legend()
    xlim(-0.2, 3)
    outputname = 'figs/walls-%02d.pdf' % (reduced_density)
    savefig(outputname, bbox_inches=0)
    print('figs/walls-%02d.pdf' % (reduced_density))


# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5]]
# savefig('figs/walls-60.pdf')
plot_walls(60, [10.0, 5.0, 2.5])

# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# savefig('figs/walls-100.pdf')
plot_walls(100, [10.0, 5.0, 2.5]) # 1.0 could work?

# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (1.5, temp) for temp in [2.5]]
# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (1.5, temp) for temp in [2.5]]
# savefig('figs/walls-150.pdf')
plot_walls(150, [2.5])
