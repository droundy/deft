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
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])

def read_walls(ff):
    filename = "mc/wallsMC-pair-0.%d0-density.dat" % (10*ff)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try a different values for ff or leave it blank to use defaults, or generate more monte carlo data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

den = read_walls(ff)

pylab.plot(den[:,0], den[:,1])
pylab.title('density(z), $ff = %g$' %(ff))
pylab.xlabel("z")
pylab.ylabel("density (balls/unit$^3$)")

pylab.show()
