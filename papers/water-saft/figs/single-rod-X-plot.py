#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

colors = ["#990022", "#dd6622", "#eecc22", "#22b544", "#3333bb"]
radii = [0.2, 0.6, 1.0]

# oldline = "--"
# newline = "-"

nm = 18.8972613
hardsphereR = 3.03420/2/10 # from Clark et al, in nm

for i in range(len(radii)):
    newdata = pylab.loadtxt('figs/single-rod-slice-%04.2f.dat' % (2*radii[i]))
    hugdata = pylab.loadtxt('figs/hughes-single-rod-slice-%04.2f.dat' % (2*radii[i]))
    rnew = newdata[:,0]/nm
    rhug = hugdata[:,0]/nm
    newcutoff = radii[i]
    hugcutoff = radii[i] - hardsphereR
    newdata[rnew<newcutoff] = 1
    pylab.plot(rnew[rnew>=hugcutoff], 4*(1-newdata[:,2])[rnew>=hugcutoff], color = colors[i], linestyle='-')
    pylab.plot(rhug[rhug>=hugcutoff], 4*(1-hugdata[:,2])[rhug>=hugcutoff], color = colors[i], linestyle='--')
    #pyplot.vlines(x = radii[i], ymin = 0, ymax = 4, color=colors[i], linestyles=':')
    pyplot.vlines(x = hugcutoff, ymin = 0, ymax = 4, color=colors[i], linestyles=':')

#plot properties
pyplot.ylabel('Number of hydrogen bonds')
pyplot.xlabel('Radius (nm)')
pyplot.xlim(0, 1.3)
pyplot.savefig('figs/single-rod-X-plot.pdf')
