#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab
import RG

# Read in data
data = pylab.loadtxt('figs/npart_RG-out.dat')

T = data[:,0]
# nvapor = data[:,1]
# nliquid = data[:,2]

etapart = data[:,5]*pylab.pi*RG.sigma**3/6

etavapor = data[:,1]*pylab.pi*RG.sigma**3/6
etaliquid = data[:,2]*pylab.pi*RG.sigma**3/6

# Plot the curve
pylab.plot(etavapor, T, 'b-')
pylab.plot(etaliquid, T, 'b-')
# pylab.plot(etapart, T, 'r-')

pylab.xlabel(r'$\eta$')
pylab.ylabel('T')
pylab.title('Liquid-Vapor Coexistence, i=1')

pylab.savefig('figs/coexistance_RG.pdf')
pylab.show()
