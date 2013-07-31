#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

newdata = pylab.loadtxt('figs/single-rod-in-water.dat')
hugdata = pylab.loadtxt('figs/hughes-single-rod-in-water.dat')
rnew = newdata[:,0]
rhug = hugdata[:,0]
newbrokenHB = newdata[:,2]
hugbrokenHB = hugdata[:,2]
pylab.plot(rhug, newbrokenHB, color = 'red', linestyle='-')
pylab.plot(rhug, hugbrokenHB, color = 'blue', linestyle='--')

#plot properties
pyplot.ylabel('Broken bonds?')
pyplot.xlabel('Radius (nm)')
pyplot.savefig('figs/single-rod-broken-HB.eps')
