#!/ur/bin/env python

import matplotlib

from matplotlib import pyplot
import pylab, numpy

ofile = numpy.loadtxt('four-rods-energy-1.0.dat')
distance = ofile[:,0]
energy = ofile[:,1]

EvsD = pyplot.figure()
pyplot.scatter(distance,energy)
pyplot.show()
