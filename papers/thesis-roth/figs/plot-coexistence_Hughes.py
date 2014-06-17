#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab
import Hughes

# Read in data
data = pylab.loadtxt('figs/npart_Hughes-out.dat')

T = data[:,0]
nvapor = data[:,1]
nliquid = data[:,2]

# Plot the curve
pylab.plot(nvapor, T,'.-')
pylab.plot(nliquid, T,'.-')

pylab.xlabel('n')
pylab.ylabel('T (K)')

pylab.savefig('figs/coexistence_Hughes.pdf')
pylab.show()
