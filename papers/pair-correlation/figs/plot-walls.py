#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
import pylab, numpy, sys
import os.path

ff = .3
z0 = 9
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])
if (len(sys.argv) > 2):
    z0 = float(sys.argv[2])



constants = numpy.loadtxt("constants.dat")
zmax = constants[0]
xmax = constants[1]
dx = constants[2]

x = numpy.arange(0, xmax, dx)
z = numpy.arange(0, zmax, dx)
Z, X = numpy.meshgrid(z, x)

def read_walls(ff, z0):
    filename = "wallsWB-pair-0.%d0-%g.dat" % (10*ff, z0)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try different values for ff and z0, or leave them blank to use defaults, or edit walls.cpp then rerun make papers to generate more data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

g2 = read_walls(ff, z0)

pylab.figure(1)
CS = pylab.contourf(Z, X, g2, 100)
CB = pylab.colorbar(CS)

pylab.axvline(x=0, color='k', linestyle=':')
pylab.axhline(y=0, color='k', linestyle=':')
pylab.title('$g^{(2)}(z_0, z_1, x_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
pylab.xlabel("z")
pylab.ylabel("x")

pylab.figure(2)
pylab.plot(z, g2[0,:])

pylab.title('$g^{(2)}(z_0, z_1, x_1)$, $z_0 = %g$, $ff = %g$ for $x_1=0$' %(z0, ff))
pylab.xlabel("z")
pylab.ylabel("g")

pylab.figure(3)
# This should be exactly the same as g(r) in contact/figs/plot-ghs2.py for the same
# gsigma, if everything is coded correctly. It is scaled similary so the two can be
# compared.
n = int(z0/dx)
pylab.plot(x, g2[:,n])

pylab.title('$g^{(2)}(z_0, z_1, x_1)$, $z_0 = %g$, $ff = %g$ for $z_1=z_0$.' %(z0, ff))
pylab.xlim(0,2.25)
pylab.ylim(0, 3.5)
pylab.xlabel("x")
pylab.ylabel("g")
pylab.show()
