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

def plot_sw(ff, temps):
    figure()
    for temp in temps:
        fname = 'data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, ff/100.0)
        data = loadtxt(fname)
        r = data[:,0]
        filling_fraction = data[:,1]
        plot(r, filling_fraction/(ff/100.0), label = 'T = %.2f' %(temp))

    title('Radial distribution function ff = %g' % (ff/100))
    xlabel(r'$r$')
    ylabel('$g(r)$')
    legend()
    xlim(-0.2, 4)
    outputname = 'figs/radial-sw-%02d.pdf' % (ff)
    savefig(outputname, bbox_inches=0)
    print('figs/radial-sw-%02d.pdf' % (ff))


# input: ['data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, 0.2) for temp in [10.0, 5.0, 3.0, 1.0]]
# savefig('figs/radial-sw-20.pdf')
plot_sw(20, [10.0, 5.0, 3.0, 1.0])

# input: ['data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, 0.2) for temp in [10.0, 5.0, 3.0, 1.0]]
# savefig('figs/radial-sw-30.pdf')
plot_sw(30, [10.0, 5.0, 3.0, 1.0])

# figure()
# X = loadtxt('radial-sw-1.00-1.30-0.30-X.dat')
# Y = loadtxt('radial-sw-1.00-1.30-0.30-Y.dat')
# gradial = loadtxt('radial-sw-1.00-1.30-0.30-n.dat')
# contourf(X,Y,gradial,100)
# colorbar()
# savefig('contour-square-well.pdf')

