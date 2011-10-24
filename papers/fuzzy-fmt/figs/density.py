#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 5:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])

pylab.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"b--",label='WB density')
pylab.plot(wbtdata[:,0],wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT density')

pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.legend(loc='upper left', ncol=2).get_frame().set_alpha(0.5)

pylab.savefig(sys.argv[4])
