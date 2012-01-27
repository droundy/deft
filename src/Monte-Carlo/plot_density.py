#!/usr/bin/python

import pylab, numpy, sys

div = 16

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)


d = div 
fignum = 0
for f in sys.argv[1:]:
    data = numpy.loadtxt(f)
    data
    pylab.figure(fignum)
    fignum+= 1
    d = len(data[:,0])/4
    pylab.plot(data[range(0,d),0],data[range(0,d),1]*4*numpy.pi/3,"o-")
    pylab.plot(data[range(0,d),0],data[range(0,d),2]*4*numpy.pi/3,"o-")#,label="ConDensity")
    pylab.plot(data[range(0,d),0],data[range(0,d),3]*4*numpy.pi/3,"o-")#,label="CenConDensity")
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend()

    pylab.figure(fignum)
    fignum+= 1
    pylab.plot(data[range(d,2*d),0],data[range(d,2*d),1]*4*numpy.pi/3,"o-")
    pylab.plot(data[range(d,2*d),0],data[range(d,2*d),2]*4*numpy.pi/3,"o-")#,label="ConDensity")
    pylab.plot(data[range(d,2*d),0],data[range(d,2*d),3]*4*numpy.pi/3,"o-")#,label="CenConDensity")
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend()

    pylab.figure(fignum)
    fignum+= 1
    pylab.plot(data[range(2*d,3*d),0],data[range(2*d,3*d),1]*4*numpy.pi/3,"o-")
    pylab.plot(data[range(2*d,3*d),0],data[range(2*d,3*d),2]*4*numpy.pi/3,"o-")#,label="ConDensity")
    pylab.plot(data[range(2*d,3*d),0],data[range(2*d,3*d),3]*4*numpy.pi/3,"o-")#,label="CenConDensity")
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend()

    pylab.figure(fignum)
    fignum+= 1
    pylab.plot(data[range(3*d,4*d),0],data[range(3*d,4*d),1]*4*numpy.pi/3,"o-")
    pylab.plot(data[range(3*d,4*d),0],data[range(3*d,4*d),2]*4*numpy.pi/3,"o-")#,label="ConDensity")
    pylab.plot(data[range(3*d,4*d),0],data[range(3*d,4*d),3]*4*numpy.pi/3,"o-")#,label="CenConDensity")
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend()

pylab.show()
