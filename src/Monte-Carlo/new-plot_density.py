#!/usr/bin/python
# -*- coding: utf-8 -*-

import pylab, numpy, sys

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)

for f in sys.argv[1:]:
    data = numpy.loadtxt(f)
    #print data[:,0]
    #print data[1]
    t= 22/19.0
    pylab.plot(data[:,0],data[:,1]*4*numpy.pi/3,"o-", label='density')
    pylab.plot(data[:,0],data[:,10]*4*numpy.pi/3,"r-", label='n0')
    pylab.plot(data[:,0],t*data[:,2]*4*numpy.pi/3,"x-", label='small contact density')
    pylab.plot(data[:,0],t*data[:,3]*4*numpy.pi/3,"o-", label='small centered contact density')
    pylab.plot(data[:,0],t*data[:,4]*4*numpy.pi/3,"x-", label='medium contact density')
    pylab.plot(data[:,0],t*data[:,5]*4*numpy.pi/3,"o-", label='medium centered contact density')
    pylab.plot(data[:,0],t*data[:,6]*4*numpy.pi/3,"x-", label='large contact density')
    pylab.plot(data[:,0],t*data[:,7]*4*numpy.pi/3,"o-", label='large centered contact density')
    pylab.plot(data[:,0],t*data[:,8]*4*numpy.pi/3,"x-", label='huge contact density')
    pylab.plot(data[:,0],t*data[:,9]*4*numpy.pi/3,"o-", label='huge centered contact density')
    pylab.xlabel("radius")
    pylab.ylabel("filling fraction")
    pylab.legend(loc='upper left', ncol=1).get_frame().set_alpha(0.5)
    #pylab.yticks( numpy.arange(0.0,5.0, 0.1) )
    #pylab.grid('True')
    pylab.show()


