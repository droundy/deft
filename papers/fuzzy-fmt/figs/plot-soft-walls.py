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
        data[i,1] = data[abs(data[:,0]) == abs(data[i,0]),1].mean()
    return data[:len(data)//2, :]

def plot_soft_walls(reduced_density, temps):
    figure()
    sigma_over_R=2**(5/6)
    have_labelled_dft = False
    NUM = 1
    for temp in temps:
        fname = 'figs/new-data/soft-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:,0]
        nreduced_density = data[:,1]
        if have_labelled_dft:
            plot(z, nreduced_density, styles.new_dft_code(temp))
        else:
            plot(z, nreduced_density, styles.new_dft_code(temp), label = 'DFT $T^* = %g$' % temp)
            have_labelled_dft = True
        
        fname = 'figs/mc-soft-wall-%04.4f-%.4f-2700.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        zmin = data[:,0].min()
        plot(smooth(data[:,0]-zmin,NUM), smooth(data[:,1],NUM)/2**(-5.0/2.0),
             styles.mcwca(temp), label = 'WCA MC $T^*$ = %g' % temp)

        fname = 'figs/new-data/bh-soft-wall-%.2f-%.2f.dat' % (reduced_density/100.0, temp)
        data = loadtxt(fname)
        z = data[:,0]
        plot(z, data[:,1], styles.bh_dft(temp))

    #plot(data[:,0], data[:,2]*0.1, 'r:', label='$V_{wall}$ (arbitrary units)')

    title('Soft walls with bulk $n^* = %g$' % (reduced_density/100))
    xlabel(r'$z/\sigma$')
    ylabel('$n^*(r)$')
    legend()
    xlim(-0.2, 6)
    outputname = 'figs/soft-walls-%02d.pdf' % (reduced_density)
    savefig(outputname, bbox_inches=0)
    print('figs/walls-%02d.pdf' % (reduced_density))


# input: ['figs/mc-soft-wall-%04.4f-%.4f-2700.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# input: ['figs/new-data/soft-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# input: ['figs/new-data/bh-soft-wall-%.2f-%.2f.dat' % (0.6, temp) for temp in [10.0, 5.0, 2.5, 1.0, 0.1]]
# savefig('figs/soft-walls-60.pdf')
plot_soft_walls(60, [10.0, 5.0, 2.5, 1.0, 0.1])
