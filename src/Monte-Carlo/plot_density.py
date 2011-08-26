#!/usr/bin/python

import pylab, numpy

data = numpy.loadtxt("density.dat")
#print data[:,0]
#print data[1]
pylab.plot(data[:,0],data[:,1]*4*numpy.pi/3,"o-")
pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.show()
