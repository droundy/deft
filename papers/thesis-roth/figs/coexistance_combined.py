#!/usr/bin/python

import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import pylab

# Read in data
data_Hughes = pylab.loadtxt('figs/npart_Hughes-out.dat')
data_SW = pylab.loadtxt('figs/npart_SW-out.dat')

T_Hughes = data_Hughes[:,0]
nvapor_Hughes = data_Hughes[:,1]
nliquid_Hughes = data_Hughes[:,2]

T_SW = data_SW[:,0]
nvapor_SW = data_SW[:,1]
nliquid_SW = data_SW[:,2]

# Plot the curve
pylab.plot(nvapor_Hughes/0.374, T_Hughes/698, 'b-')
pylab.plot(nliquid_Hughes/0.374, T_Hughes/698, 'b-')

pylab.plot(nvapor_SW/0.0737, T_SW/8.6, 'r-')
pylab.plot(nliquid_SW/0.0737, T_SW/8.6, 'r-')

pylab.xlabel(r'$n/n_{crit}$')
pylab.ylabel(r'$T/T_{crit}$')

pylab.savefig('figs/combined-coexistence-plot.pdf')
pylab.show()
