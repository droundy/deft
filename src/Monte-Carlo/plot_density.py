#!/usr/bin/python

import pylab, numpy

data = numpy.loadtxt("density.dat")
#print data[:,0]
#print data[1]
pylab.plot(data[:,0],data[:,1]*4*numpy.pi/3,"o-")
pylab.plot(data[:,0],data[:,2]*4*numpy.pi/3,"o-",label="ConDensity")
pylab.plot(data[:,0],data[:,3]*4*numpy.pi/3,"o-",label="CenConDensity")
pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.legend()
pylab.show()
