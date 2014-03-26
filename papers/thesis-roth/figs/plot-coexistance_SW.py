#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab
import SW

# Read in data
data = pylab.loadtxt('figs/npart_SW-out.dat')

T = data[:,0]
nvapor = data[:,1]
nliquid = data[:,2]

etapart = data[:,5]*pylab.pi*SW.sigma**3/6

etavapor = data[:,1]*pylab.pi*SW.sigma**3/6
etaliquid = data[:,2]*pylab.pi*SW.sigma**3/6

# Plot the curve
pylab.plot(etavapor, T, 'b-')
pylab.plot(etaliquid, T, 'b-')
# pylab.plot(etapart, T, 'r-')

pylab.xlabel(r'$\eta$')
pylab.ylabel('T')

pylab.savefig('figs/coexistance_SW.pdf')
pylab.show()


