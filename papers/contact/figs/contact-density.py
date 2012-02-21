#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

if len(sys.argv) != 6:
    print("Usage:  " + sys.argv[0] + " mc-filename.dat wb-filename.dat wbt-filename.dat wb-m2.dat out-filename.pdf")
    exit(1)

mcdata = numpy.loadtxt(sys.argv[1])
dftdata = numpy.loadtxt(sys.argv[2])
wbtdata = numpy.loadtxt(sys.argv[3])
wbm2data = numpy.loadtxt(sys.argv[4])

pylab.plot(mcdata[:,0],mcdata[:,1]*4*numpy.pi/3,"b-",label='MC density')
pylab.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"b--",label='DFT density')
pylab.plot(wbtdata[:,0],wbtdata[:,1]*4*numpy.pi/3,"m-.",label='WBT density')
pylab.plot(wbm2data[:,0],wbm2data[:,1]*4*numpy.pi/3,"c--",label='WB mark II density')

pylab.plot(dftdata[:,0],dftdata[:,6]*4*numpy.pi/3,"c--",label="$n_0$")

off = 2
pylab.plot(mcdata[:,0],mcdata[:,2+2*off]*4*numpy.pi/3,"r-",label="ConDensity")
pylab.plot(mcdata[:,0],mcdata[:,3+2*off]*4*numpy.pi/3,"g-",label="CenConDensity")

me = 3
#pylab.plot(dftdata[:,0],dftdata[:,3]*4*numpy.pi/3,"g+--",label="simple contact")
pylab.plot(dftdata[:,0],dftdata[:,4]*4*numpy.pi/3,"gx--",label="Yu and Wu")
pylab.plot(dftdata[:,0],dftdata[:,5]*4*numpy.pi/3,"ro--",label="DFT at sphere")
#pylab.plot(wbm2data[:,0],wbm2data[:,5]*4*numpy.pi/3,"r+--",markevery=me,label="DFT at sphere (mark II)",
#           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)
pylab.plot(dftdata[:,0],dftdata[:,7]*4*numpy.pi/3,"rx--",markevery=me,label="Gross",
           markerfacecolor='none',markeredgecolor='red', markeredgewidth=1)

pylab.xlabel("radius")
pylab.ylabel("filling fraction")
pylab.legend(loc='upper left', ncol=2).get_frame().set_alpha(0.5)

pylab.savefig(sys.argv[5])
