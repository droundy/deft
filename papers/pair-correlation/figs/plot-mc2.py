#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
#matplotlib.use('Agg')
import pylab, numpy, sys
import os.path
import matplotlib.colors as mcolors

ff = .3
z0 = 2
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])
if (len(sys.argv) > 2):
    z0 = float(sys.argv[2])

zmax = 10
rmax = 5

def read_walls(ff, z0):
    filename = "mc2/wallsMC-pair-0.%d0-%1.2f.dat" % (10*ff, z0)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try different values for ff and z0, or leave them blank to use defaults, or generate more monte carlo data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

g2 = read_walls(ff, z0)

#g2*= 20.25

zbins = len(g2[0,:])
rbins = len(g2[:,0])
dr = rmax/rbins
dz = zmax/zbins

r = numpy.arange(0, rmax, dr)
z = numpy.arange(0, zmax, dz)
Z, R = numpy.meshgrid(z, r)

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
map = mcolors.LinearSegmentedColormap('map', cdict)
pylab.register_cmap(cmap = map)
cmap = pylab.get_cmap('map')

fig = pylab.figure(1)
max = int(g2.max()) + 1
print max
dw = 0.05
levels = numpy.arange(2-max, max, dw)

#CS = pylab.contourf(Z, R, g2, levels, cmap=cmap)
#pylab.contourf(Z, -R, g2, levels, cmap=cmap)


CS = pylab.contourf(Z, R, g2, 50)
pylab.contourf(Z, -R, g2, 50)

CB = pylab.colorbar(CS)
pylab.axes().set_aspect('equal')

pylab.title('$g^{(2)}(z_0, z_1, r_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
pylab.xlabel("z")
pylab.ylabel("r")

pylab.show()
