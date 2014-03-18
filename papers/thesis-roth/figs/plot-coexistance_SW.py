#!/usr/bin/python

import pylab

# Read in data
data = pylab.loadtxt('npart_SW-out.dat')

T = data[:,0]
nvapor = data[:,1]
nliquid = data[:,2]

etavapor = data[:,1]*pylab.pi*2**3/6
etaliquid = data[:,2]*pylab.pi*2**3/6

# Plot the curve
pylab.plot(etavapor, T)
pylab.plot(etaliquid, T)

pylab.xlabel(r'$\eta$')
pylab.ylabel('T')

pylab.show()
