#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab
import SW

# Read in data
data = pylab.loadtxt('npart_SW-out.dat')

T = data[:,0]
nvapor = data[:,1]
nliquid = data[:,2]

etavapor = data[:,1]*pylab.pi*SW.sigma**3/6
etaliquid = data[:,2]*pylab.pi*SW.sigma**3/6

# Plot the curve
pylab.plot(etavapor, T)
pylab.plot(etaliquid, T)

pylab.xlabel(r'$\eta$')
pylab.ylabel('T')

pylab.savefig('coexistance_SW.pdf')
# pylab.show()


