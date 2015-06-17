#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if not 'show' in sys.argv:
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
    #ylim(0, 4)
    outputname = 'figs/radial-sw-%02d.pdf' % (ff)
    savefig(outputname, bbox_inches=0)
    print('figs/radial-sw-%02d.pdf' % (ff))


# input: ['data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, 0.2) for temp in [10.0, 5.0, 3.0]]
# savefig('figs/radial-sw-20.pdf')
plot_sw(20, [10.0, 5.0, 3.0])

# input: ['data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, 0.3) for temp in [10.0, 5.0, 3.0]]
# savefig('figs/radial-sw-30.pdf')
plot_sw(30, [10.0, 5.0, 3.0])

figure()
#col=matplotlib.cm.get_cmap('spring')
#col.set_under('w')
X = loadtxt('data/radial-sw-3.00-1.30-0.20-X.dat')
Y = loadtxt('data/radial-sw-3.00-1.30-0.20-Y.dat')
gradial = loadtxt('data/radial-sw-3.00-1.30-0.20-eta.dat')
gradial /= 0.20

gmax = gradial.max()
xlo = 0.5/gmax
xwhite = 1.0/gmax
xhi = 0.3*xwhite + 0.7
xhier = (1 + xhi)/2.0
print [xlo, xwhite, xhi, xhier]

cdict = {'red':   [(0.0,  0.0, 0.0),
                   (xlo,  1.0, 1.0),
                   (xwhite,  1.0, 1.0),
                   (xhi,  0.0, 0.0),
                   (xhier,0.0, 0.0),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0, 0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (xwhite, 1.0, 1.0),
                   (xhi, 0.0, 0.0),
                   (xhier,1.0, 1.0),
                   (1.0, 1.0, 1.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (xwhite,  1.0, 1.0),
                   (xhi,  1.0, 1.0),
                   (xhier,0.0, 0.0),
                   (1.0,  0.0, 0.0)]}
cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)

pcolormesh(X,Y,gradial, vmax=gmax, vmin=0, cmap=cmap)
axis('tight')
axes().set_aspect('equal')
colorbar()
savefig('figs/contour-square-well.pdf')

show()
