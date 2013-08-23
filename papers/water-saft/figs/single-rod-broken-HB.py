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
p1, = pylab.plot(rnew, newbrokenHB/2, color = '#990022', linestyle = '-')  #'r-')
p2, = pylab.plot(rhug, hugbrokenHB/2, color = '#220099', linestyle = '--') #'k--')

#plot properties
pyplot.ylabel('Broken bonds per nm')
pyplot.xlabel('Radius of rod (nm)')
pyplot.xlim(0, 1.5)
pyplot.ylim(0, 10)
pyplot.legend([p2,p1],["Hughes, et al", "This work"], loc = 2)
pyplot.savefig('figs/single-rod-broken-HB.pdf')
