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
import os, readandcompute

colordict = { 3.0: 'b', 5.0: 'g', 10.0: 'r' }
def color(T):
    if T in colordict.keys():
        return colordict[T]
    return 'k'

def plot_sw(ff, temps):
    figure()
    for temp in temps:
        fname = 'data/radial-sw-%.2f-%.2f-%.2f.dat' % (temp, 1.3, ff/100.0)
        data = loadtxt(fname)
        r = data[:,0]
        filling_fraction = data[:,1]
        plot(r, filling_fraction/(ff/100.0), color(temp)+'-',
             label = 'T = %.2f' %(temp))

        N = 500
        ww = 1.3
        g, r = readandcompute.g_r('data/mc/ww%.2f-ff%.2f-N%d' % (ww,ff/100.0,N), temp)
        plt.plot(r/2, g, color(temp)+':')

    title('Radial distribution function ff = %g' % (ff/100))
    xlabel(r'$r$')
    ylabel('$g(r)$')
    legend()
    xlim(-0.2, 4)
    #ylim(0, 4)

plot_sw(20, [10.0, 5.0, 3.0])
plt.savefig('figs/radial-sw-20.pdf')

plot_sw(30, [10.0, 5.0, 3.0])
plt.savefig('figs/radial-sw-30.pdf')

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
