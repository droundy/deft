#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
import pylab, numpy, sys
import os.path
from matplotlib.colors import LinearSegmentedColormap

ff = .3
z0 = 2
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])
if (len(sys.argv) > 2):
    z0 = float(sys.argv[2])



constants = numpy.loadtxt("constants.dat")
zmax = 20
xmax = 10
bins = 200

dx = xmax/bins
dz = zmax/bins

x = numpy.arange(0, xmax, dx)
z = numpy.arange(0, zmax, dz)
Z, X = numpy.meshgrid(z, x)

def read_walls(ff, z0):
    filename = "mc/wallsMC-pair-0.%d0-%1.1f.dat" % (10*ff, z0)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try different values for ff and z0, or leave them blank to use defaults, or generate more monte carlo data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

g2 = read_walls(ff, z0)
#g2 /= 1.25
#g2*=.76
cdict = {'red':  ((0.0, 0.0, 0.0),
                  (0.25,0.0, 0.0),
                  (0.5, 0.8, 1.0),
                  (0.75,1.0, 1.0),
                  (1.0, 0.4, 1.0)),

        'green': ((0.0, 0.0, 0.0),
                  (0.25,0.0, 0.0),
                  (0.5, 0.9, 0.9),
                  (0.75,0.0, 0.0),
                  (1.0, 0.0, 0.0)),

        'blue':  ((0.0, 0.0, 0.4),
                  (0.25,1.0, 1.0),
                  (0.5, 1.0, 0.8),
                  (0.75,0.0, 0.0),
                  (1.0, 0.0, 0.0))
        }
map = LinearSegmentedColormap('map', cdict)
pylab.register_cmap(cmap = map)
cmap = pylab.get_cmap('map')
fig = pylab.figure(1)
max = int(g2.max()) + 1
print max
dw = .05
levels = numpy.arange(2-max, max+dw/2, dw)

CS = pylab.contourf(Z, X, g2)#, levels, cmap=cmap)
pylab.contourf(Z, -X, g2)#, levels, cmap=cmap)
CB = pylab.colorbar(CS)
pylab.axes().set_aspect('equal')

pylab.title('$g^{(2)}(z_0, z_1, x_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
pylab.xlabel("z")
pylab.ylabel("x")


pylab.show()
