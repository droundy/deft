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
    NUM=1
    for temp in temps:
        fname = 'figs/new-data/wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:, 0]
        z -= z.min() + 3
        nreduced_density = data[:, 1]
        if have_labelled_dft:
            plot(z, nreduced_density, styles.new_dft_linestyle(), color = styles.color_from_kT(temp))
        else:
            plot(z, nreduced_density, styles.new_dft_linestyle(), color = styles.color_from_kT(temp),
                 label = 'DFT $T^* = %g$' % temp)
            have_labelled_dft = True

        fname = 'figs/new-data/bh-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:, 0]
        z -= z.min() + 3
        nreduced_density = data[:, 1]
        if have_labelled_bh:
            plot(z, nreduced_density, styles.bh_dft_linestyle(), color = styles.color[temp])
        else:
            plot(z, nreduced_density, styles.bh_dft_linestyle(), color = styles.color[temp], label = 'BH $T^* = %g$' % temp)
            have_labelled_bh = True

        fname = 'figs/mcfcc-walls-%04.4f-%.4f.dat' % (reduced_density/100.0, temp)
        if os.path.exists(fname):
            data = loadtxt(fname)
            z = data[:,0]
            z -= z.min()
            smoothed = smooth(data[:, 1], NUM)
            print(reduced_density/100, temp, smoothed[len(smoothed)//2])
            plot(smooth(z, NUM)/sigma_over_R, smooth(data[:, 1], NUM)*sigma_over_R**3,
                styles.mcwca_linestyle(), color = styles.color_from_kT(temp), label = 'WCA MC $T^*$ = %g' % temp)

    title('Hard walls with bulk $n^* = %g$' % (reduced_density/100), fontsize=20)
    xlabel(r'$z/\sigma$', fontsize=18)
    plt.xticks(fontsize=17)
    ylabel('$n^*(r)$', fontsize=18)
    plt.yticks(fontsize=17)
    legend(fontsize=18)
    plt.ylim(0)
    plt.xlim(-.25, 5)

# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# input: ['figs/new-data/bh-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0]]

# input: ['figs/new-data/wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/bh-wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]

# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (0.6, temp) for temp in [10.0, 2.5]]
# input: ['figs/mcfcc-walls-%04.4f-%.4f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]

subplot(2, 1, 1)
plot_walls(60, [10.0, 2.5, 1.0])

subplot(2, 1, 2)
plot_walls(100, [10.0, 5.0, 2.5]) # 1.0 could work?

plt.tight_layout();
savefig('figs/hard-walls.pdf', bbox_inches=0)
print('figs/hard-walls.pdf')

