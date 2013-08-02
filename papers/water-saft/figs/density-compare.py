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
blueish = '#55dae0'

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

hugdata = pylab.loadtxt('figs/hughes-single-rod-slice-0.40.dat')
rhug = hugdata[:,0]/nm
hugdensity = hugdata[:,1]/gpermL
pylab.plot(rhug, hugdensity, color = '#2222aa', linestyle='--')

newdata = pylab.loadtxt('figs/single-rod-slice-0.40.dat')
rnew = newdata[:,0]/nm
newdensity = newdata[:,1]/gpermL
pylab.plot(rnew, newdensity, color = '#118844', linestyle='-')

pyplot.hlines(1, 0, 1.3, 'black', ':')

circleheight = 0.25
ymax = 3.9
rmax = 1.2
hardsphere_diameter = 3.0342/10 # nm
rod_radius = 0.2 # nm

xpoints = [rod_radius + n*hardsphere_diameter for n in range(4)]
ypoints = [circleheight]*4
pyplot.plot(xpoints, ypoints, marker = 'o', color = 'black', linestyle = '')
fig = pyplot.gcf()
for n in range(4):
    xpos = rod_radius + n*hardsphere_diameter
    pyplot.vlines(xpos, 0, ymax, grey, ':')
    fig.gca().add_artist(Ellipse((xpos, circleheight),
                                 hardsphere_diameter, 1.2*hardsphere_diameter*ymax/rmax,
                                 color = blueish, fill=False))

#plot properties
pyplot.ylabel('Density (g/mL)')
pyplot.xlabel('Radius (nm)')
pyplot.ylim(0, ymax)
pyplot.xlim(0, rmax)
pyplot.savefig('figs/density-compare.pdf')
