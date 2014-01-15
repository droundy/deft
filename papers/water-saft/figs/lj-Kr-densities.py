#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab
from matplotlib.patches import Ellipse

import Bowron_Kr_320K
import Bowron_Kr_278K

nm = 18.8972613
gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density
rmax = 1.3

grey = '#999999'
blueish = '#99cccc'

temps = [  278, 320 ]
for t in temps:
    pylab.figure()
    data = pylab.loadtxt('figs/lj-Kr-%gK-new-precond.dat' %t)
    data_hughes = pylab.loadtxt('figs/hughes-lj-Kr-%gK.dat' %t)
    r = data[:,1]/nm
    r_hughes = data_hughes[:,1]/nm
    density = data[:,3]/gpermL
    dens_hughes = data_hughes[:,3]/gpermL
    pylab.plot(r_hughes*10, dens_hughes, color = grey, linestyle='--', label='Hughes, et al')
    pylab.plot(r*10, density, color = '#3333aa', linestyle='-', label='This work')

    if t == 320:
        pylab.plot(Bowron_Kr_320K.r, Bowron_Kr_320K.g, 'k:', label='experiment')
    if t == 278:
        pylab.plot(Bowron_Kr_278K.r, Bowron_Kr_278K.g, 'k:', label='experiment')
    pylab.title('Kr at %g K and 100 bars' % t)

    pyplot.ylabel('$g_{O-Kr}$')
    pyplot.xlabel(r'Radius ($\AA$)')
    pyplot.legend()
    pyplot.xlim(2, 8)
    if t == 320:
        pylab.savefig('figs/Kr-320K-densities.pdf')
    if t == 278:
        pylab.savefig('figs/Kr-278K-densities.pdf')

pyplot.show()
