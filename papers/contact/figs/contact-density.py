#!/usr/bin/python

import pylab, numpy, sys

if len(sys.argv) != 4:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat dft-filename.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])

pylab.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT density')

pylab.plot(mcdata[:,0],mcdata[:,2]*4*numpy.pi/3,"r-",label="ConDensity")
pylab.plot(mcdata[:,0],mcdata[:,3]*4*numpy.pi/3,"g-",label="CenConDensity")

pylab.plot(dftdata[:,0],dftdata[:,3]*4*numpy.pi/3,"g+--",label="simple contact")
pylab.plot(dftdata[:,0],dftdata[:,4]*4*numpy.pi/3,"gx--",label="Wu and Fu")
pylab.plot(dftdata[:,0],dftdata[:,5]*4*numpy.pi/3,"ro--",label="DFT at sphere")

pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.legend(loc='upper left', ncol=2)

pylab.savefig(sys.argv[3])
