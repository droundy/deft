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

zmax = 20
rmax = 10

def read_walls(ff, z0):
    filename = "mc/wallsMC-pair-0.%d-%1.2f.dat" % (10*ff, z0)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try different values for ff and z0, or leave them blank to use defaults, or generate more monte carlo data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

g2 = read_walls(ff, z0)


rbins = len(g2)
dr = rmax/rbins

r = numpy.arange(0, rmax, dr)
fig = pylab.figure(1)
max = int(g2.max()) + 1
print max
dw = 0.05

pylab.plot(r, g2)

pylab.title('$g^{(2)}(z_0, z_1, r_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
pylab.xlabel("r")
pylab.ylabel("g")

plotname = "radial-" + str(10*ff) + "-" + str(z0) + ".pdf"
pylab.savefig(plotname)
