#!/usr/bin/env python

#need this to run without xserver
import matplotlib
#matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

newcalc = ""
hughescalc = "hughes-"

nm = 18.8972613
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

colors = ["#44dd55", "#3377aa", "#002288"]
radii = [ 0.2, 0.6, 1.0 ]

for i in range(len(radii)):
    newdata = pylab.loadtxt('figs/single-rod-slice-%04.2f.dat' % (2*radii[i]))
    hugdata = pylab.loadtxt('figs/hughes-single-rod-slice-%04.2f.dat' % (2*radii[i]))
    rnew = newdata[:,0]/nm
    rhug = hugdata[:,0]/nm
    newdensity = newdata[:,1]/gpermL
    newHB = 4*(1-newdata[:,2])
    #hugHB = 4*(1-hugdata[:,2])
    bulkHB = 4*(1-newdata[len(newdata) - 2,2])
    #brokenHB = (bulkHB - newHB)*newdensity
    #brokenHB[newdensity<0.01] = 0
    newdenstotal = 0
    #for j in range(len(rnew)):
    #    newdenstotal += math.pi*(rnew[j+1]**2 - rnew[j]**2)*newdensity
    pylab.plot(rnew, brokenHB, color = colors[i], linestyle='-')
    #pylab.plot(rnew, newdensity, color = colors[i], linestyle='--')
    #pylab.plot(rhug, hugHB*hugdensity, color = colors[i], linestyle='--')

#plot properties
pyplot.ylabel('')
pyplot.xlabel('Radius (nm)')
pyplot.xlim(0, 1.3)
pyplot.savefig('figs/single-rod-HB-density.pdf')
pyplot.show()
