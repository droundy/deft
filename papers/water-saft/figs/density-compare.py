#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab
from matplotlib.patches import Ellipse

nm = 18.8972613
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density

grey = '#999999'
blueish = '#11aadd'

#colors = ["#44dd55", "#3377aa", "#002288"]
#radii = [ 0.2, 0.6, 1.0 ]

# for i in range(len(radii)):
#     newdata = pylab.loadtxt('figs/single-rod-slice-%04.2f.dat' % (2*radii[i]))
#     hugdata = pylab.loadtxt('figs/hughes-single-rod-slice-%04.2f.dat' % (2*radii[i]))
#     rnew = newdata[:,0]/nm
#     rhug = hugdata[:,0]/nm
#     newdensity = newdata[:,1]/gpermL
#     hugdensity = hugdata[:,1]/gpermL
#     pylab.plot(rnew, newdensity, color = colors[i], linestyle='-')
#     pylab.plot(rhug, hugdensity, color = colors[i], linestyle='--')

newdata = pylab.loadtxt('figs/single-rod-slice-0.40.dat')
rnew = newdata[:,0]/nm
newdensity = newdata[:,1]/gpermL
pylab.plot(rnew, newdensity, color = '#22bb66', linestyle='-')

pyplot.hlines(1, 0, 1.3, 'black', '--')
pyplot.vlines(0.2, 0, 4.6, grey, '--')
pyplot.vlines(0.49, 0, 4.6, grey, '--')
pyplot.vlines(0.78, 0, 4.6, grey, '--')
pyplot.vlines(1.07, 0, 4.6, grey, '--')

xpoints = [0.2, 0.49, 0.78, 1.07]
ypoints = [4.1, 4.1, 4.1, 4.1]
pyplot.plot(xpoints, ypoints, marker = 'o', color = 'black', linestyle = '')
fig = pyplot.gcf()
fig.gca().add_artist(Ellipse((0.2, 4.1), 0.29, 1.3, color = blueish, fill=False))
fig.gca().add_artist(Ellipse((0.49,4.1), 0.29, 1.3, color = blueish, fill=False))
fig.gca().add_artist(Ellipse((0.78,4.1), 0.29, 1.3, color = blueish, fill=False))
fig.gca().add_artist(Ellipse((1.07,4.1), 0.29, 1.3, color = blueish, fill=False))
#plot properties
pyplot.ylabel('Density (g/mL)')
pyplot.xlabel('Radius (nm)')
pyplot.ylim(0, 4.8)
pyplot.xlim(0, 1.3)
pyplot.savefig('figs/density-compare.pdf')
