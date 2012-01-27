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

#pylab.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0]-2,dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT density')
pylab.plot(wbtdata[:,0]-2,wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT density')

#pylab.plot(mcdata[:,0],mcdata[:,2]*4*numpy.pi/3,"r-",label="ConDensity")
#pylab.plot(mcdata[:,0],mcdata[:,3]*4*numpy.pi/3,"g-",label="CenConDensity")

me = 20
pylab.plot(dftdata[:,0]-2,dftdata[:,3]*4*numpy.pi/3,"g+--",markevery=me,label="simple contact", markeredgewidth=1)
pylab.plot(dftdata[:,0]-2,dftdata[:,4]*4*numpy.pi/3,"gx--",markevery=me,label="Yu and Wu", markeredgewidth=1)
pylab.plot(dftdata[:,0]-2,dftdata[:,5]*4*numpy.pi/3,"ro--",markevery=me,label="DFT at sphere",markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
pylab.plot(dftdata[:,0]-2,dftdata[:,6]*4*numpy.pi/3,"g*--",markevery=me,label="n2-only")

pylab.xlabel("z")
pylab.ylabel("filling fraction")
pylab.legend(loc='lower right', ncol=2).get_frame().set_alpha(0.5)
pylab.xlim(-1,8)

pylab.savefig(sys.argv[4])
