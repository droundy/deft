#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

hardsphereR = 3.03420/2/10 # from Clark et al, in nm

newdata = pylab.loadtxt('figs/sphere.dat')
hugdata = pylab.loadtxt('figs/hughes-sphere.dat')
rnew = newdata[:,0] - hardsphereR
rhug = hugdata[:,0] - hardsphereR
newbrokenHB = newdata[:,5]
hugbrokenHB = hugdata[:,5]
pylab.plot(rnew, newbrokenHB/2, color = '#990022', linestyle = '-')  #'r-')
pylab.plot(rhug, hugbrokenHB/2, color = '#220099', linestyle = '--')  #'r--')

#plot properties
pyplot.ylabel('Numer of broken bonds')
pyplot.xlabel('Radius of sphere(nm)')
pyplot.xlim(0, 1.5)
pyplot.ylim(0, 4)
pyplot.savefig('figs/sphere-broken-HB.pdf')
pyplot.show()
