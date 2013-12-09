#!/usr/bin/env python

#need this to run without xserver
import matplotlib
#matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab
from matplotlib.patches import Ellipse

nm = 18.8972613
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density
rmax = 1.3

grey = '#999999'
blueish = '#99cccc'

atoms = [ 'Xe']#'Ne', 'Ar', 'Kr', 'Xe' ]
temps = [  300 ]# 325, 350 ]
for a in atoms:
    for t in temps:
        data = pylab.loadtxt('figs/lj-%s-%gK.dat' % (a, t))
        data_hires = pylab.loadtxt('figs/lj-Xe-300K-hires.dat')
        r = data[:,1]/nm
        r_hires = data_hires[:,1]/nm
        density = data[:,3]/gpermL
        dens_hires = data_hires[:,3]/gpermL
        pylab.plot(r, density, color = '#3333aa', linestyle='-')
        pylab.plot(r_hires, dens_hires,color = grey, linestyle='--')

        #plot properties
        pyplot.ylabel('Density (g/mL)')
        pyplot.xlabel('Radius (nm)')
        #pyplot.ylim(0, ymax)
        pyplot.xlim(0, rmax)
        #pyplot.savefig('figs/lj-%s-%gK-density.pdf' % (a, t))
        pyplot.show()

        pyplot.clf()
