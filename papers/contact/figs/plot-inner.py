#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 5:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat out-filename.pdf")
    print "got len(sys.argv) of ", len(sys.argv)
    exit(1)

rmin = 0
mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])

#pylab.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0]-rmin,dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT density')
pylab.plot(wbtdata[:,0]-rmin,wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT density')

#pylab.plot(mcdata[:,0],mcdata[:,2]*4*numpy.pi/3,"r-",label="ConDensity")
#pylab.plot(mcdata[:,0],mcdata[:,3]*4*numpy.pi/3,"g-",label="CenConDensity")

pylab.plot(dftdata[:,0]-rmin,dftdata[:,3]*4*numpy.pi/3,"g+--",label="simple contact")
pylab.plot(dftdata[:,0]-rmin,dftdata[:,4]*4*numpy.pi/3,"gx--",label="Yu and Wu")
pylab.plot(dftdata[:,0]-rmin,dftdata[:,5]*4*numpy.pi/3,"ro--",label="DFT at sphere")
pylab.plot(dftdata[:,0]-rmin,dftdata[:,6]*4*numpy.pi/3,"g*--",label="n2-only")

pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.legend(loc='lower right', ncol=2).get_frame().set_alpha(0.5)
pylab.xlim(0,8)

pylab.savefig(sys.argv[4])
