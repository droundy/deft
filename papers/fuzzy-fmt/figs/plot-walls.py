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

def smooth(x, N):
    '''
    smooth(x,N) takes a 2D array x that has many columns, and averages
    out N nearby points.
    '''
    n0 = len(x)
    n = n0 - n0 % N
    x = x[:n]
    y = zeros_like(x[0::N])
    for i in range(N):
        y += x[i::N]
    return y/N

def plot_walls(reduced_density, temps):
    num = 1
    sigma_over_R=2**(5/6)
    have_labelled_dft = False
    have_labelled_bh = False
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

        fname = 'figs/new-data/bh-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:,0]
        nreduced_density = data[:,1]
        if have_labelled_bh:
            plot(z, nreduced_density, styles.color[temp]+':')
        else:
            plot(z, nreduced_density, styles.color[temp]+':', label = 'BH $T^* = %g$' % temp)
            have_labelled_bh = True
        
        fname = 'figs/mcfcc-walls-%04.4f-%.4f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        zmin = data[:,0].min()
        plot(smooth((data[:,0]-zmin)/sigma_over_R,num), smooth(data[:,1]*sigma_over_R**3,num),
             styles.mcwca(temp), label = 'WCA MC $T^*$ = %g' % temp)

    title('Hard walls with bulk $n^* = %g$' % (reduced_density/100))
    xlabel(r'$z/\sigma$')
    ylabel('$n^*(r)$')
    legend()
    xlim(-0.2, 4)

# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (0.6, temp) for temp in [10.0, 2.5, 0.1]]
# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 2.5, 0.1]]
# input: ['figs/new-data/bh-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 2.5, 0.1]]
# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/bh-wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# savefig('figs/hard-walls.pdf')
figure(figsize=(6,11))

subplot(2, 1, 1)
plot_walls(60, [10.0, 2.5, 0.1])

subplot(2, 1, 2)
plot_walls(100, [10.0, 5.0, 2.5]) # 1.0 could work?

outputname = 'figs/hard-walls.pdf'
savefig(outputname, bbox_inches=0)
print('figs/hard-walls.pdf')
