#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from scipy.special import erf
import matplotlib.lines as mlines
import os
import numpy

import styles

def smooth(x, N):
    '''
    smooth(x,N) takes a 2D array x that has many columns, and averages
    out N nearby points.
    '''
    n0 = len(x)
    n = n0 - n0 % N
    x = x[:n]
    y = numpy.zeros_like(x[0::N])
    for i in range(N):
        y += x[i::N]
    return y/N

def average_positive_and_negative(data):
    for i in range(len(data)//2):
        data[i, 1] = data[abs(data[:, 0]) == abs(data[i, 0]), 1].mean()
    return data[:len(data)//2,:]

def find_z_with_V(z, Vext, V_desired):
    z_center_index = np.abs(Vext - V_desired).argmin()
    return np.interp(V_desired, Vext[z_center_index-2:z_center_index+2], z[z_center_index-2:z_center_index+2])

def plot_soft_walls(reduced_density, temps):
    have_labelled_bh = False
    sigma_over_R=2**(5/6)
    have_labelled_dft = False
    NUM = 1
    for temp in temps:
        fname = 'figs/new-data/soft-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:, 0]
        Vext = data[:, 2]
        z_center = find_z_with_V(z, Vext, 10)
        nreduced_density = data[:, 1]
        Vext = data[:, 2]
        if have_labelled_dft:
            plot(z-z_center, nreduced_density, styles.new_dft_code(temp))
            # plot(z, Vext, 'r-')
        else:
            plot(z-z_center, nreduced_density, styles.new_dft_code(temp), label = 'DFT $T^* = %g$' % temp)
            have_labelled_dft = True

        fname = 'figs/new-data/bh-soft-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:, 0]
        Vext = data[:, 2]
        z_center = find_z_with_V(z, Vext, 10)
        nreduced_density = data[:, 1]
        if have_labelled_bh:
            plot(z-z_center, nreduced_density, styles.color[temp]+':')
        else:
            plot(z-z_center, nreduced_density, styles.color[temp]+':', label = 'BH $T^* = %g$' % temp)
            have_labelled_bh = True
        # plot(z-z_center, Vext, '--', label=f'Vext bh {temp}')


        fname = 'figs/mc-soft-wall-%04.4f-%.4f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:,0]
        zmin = z.min()
        smoothed = smooth(data[:, 1], NUM)
        Vext = data[:, 2]
        z_center = find_z_with_V(z, Vext, 10)
        print(reduced_density/100, temp, smoothed[len(smoothed)//2])
        plot(smooth(z-z_center, NUM), smooth(data[:, 1], NUM),
             styles.mcwca(temp), label = 'WCA MC $T^*$ = %g' % temp)
        # plot(z-z_center, Vext, '.', label=f'Vext mc {temp}')
        plt.ylim(0,5)

    #plot(data[:,0], data[:,2]*0.1, 'r:', label='$V_{wall}$ (arbitrary units)')

    title('Soft walls with bulk $n^* = %g$' % (reduced_density/100))
    xlabel(r'$z/\sigma$')
    ylabel('$n^*(r)$')
    legend()
    xlim(-0.3, 5)
    outputname = 'figs/soft-walls-%02d.pdf' % (reduced_density)
    savefig(outputname, bbox_inches=0)
    print(('figs/walls-%02d.pdf' % (reduced_density)))


# input: ['figs/mc-soft-wall-%04.4f-%.4f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# input: ['figs/new-data/soft-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# input: ['figs/new-data/bh-soft-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0]]
# input: ['figs/mc-soft-wall-%04.4f-%.4f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/soft-wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
# input: ['figs/new-data/bh-soft-wall-%.2f-%.2f.dat' % (1.0, temp) for temp in [10.0, 5.0, 2.5]]
figure(figsize=(9, 14.5))
subplot(2, 1, 1)
plot_soft_walls(60, [10.0, 5.0, 2.5, 1.0])

subplot(2, 1, 2)
plot_soft_walls(100, [10.0, 5.0, 2.5])

plt.tight_layout();
savefig('figs/soft-walls.pdf', bbox_inches=0)
print('figs/soft-walls.pdf')
