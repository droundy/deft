#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

hardsphereR = 3.03420/2/10 # from Clark et al, in nm

newdata = pylab.loadtxt('figs/single-rod-in-water.dat')
hugdata = pylab.loadtxt('figs/hughes-single-rod-in-water.dat')
rnew = newdata[:,0] - hardsphereR
rhug = hugdata[:,0] - hardsphereR
newbrokenHB = newdata[:,2]
hugbrokenHB = hugdata[:,2]
pylab.plot(rhug, newbrokenHB, color = 'red', linestyle='-')
pylab.plot(rhug, hugbrokenHB, color = 'blue', linestyle='--')

#plot properties
pyplot.ylabel('Broken bonds?')
pyplot.xlabel('Radius (nm)')
pyplot.xlim(0, 1.6)
pyplot.savefig('figs/single-rod-broken-HB.eps')
