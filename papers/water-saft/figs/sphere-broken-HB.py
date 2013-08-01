#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

#newdata = pylab.loadtxt('figs/sphere.dat')
hugdata = pylab.loadtxt('figs/hughes-sphere.dat')
#rnew = newdata[:,0]
rhug = hugdata[:,0]
#newbrokenHB = newdata[:,5]
hugbrokenHB = hugdata[:,5]
#pylab.plot(rhug, newbrokenHB, color = 'red', linestyle='-')
pylab.plot(rhug, hugbrokenHB, color = 'blue', linestyle='--')

#plot properties
pyplot.ylabel('Broken bonds?')
pyplot.xlabel('Radius (nm)')
pyplot.savefig('figs/sphere-broken-HB.pdf')
