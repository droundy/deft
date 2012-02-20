#!/usr/bin/python

import pylab, numpy, sys

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)

for f in sys.argv[1:]:
    data = numpy.loadtxt(f)
    #print data[:,0]
    #print data[1]
    pylab.plot(data[:,0],data[:,1]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,2]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,3]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,4]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,5]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,6]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,7]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,8]*4*numpy.pi/3,"o-")
    pylab.plot(data[:,0],data[:,9]*4*numpy.pi/3,"o-")
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend()
    pylab.show()


