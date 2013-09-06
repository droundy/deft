#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from __future__ import division
import matplotlib
matplotlib.use('Agg')

import pylab, numpy, sys

mNpermeter = 6.4230498e-07 # in atomic units
nm = 18.8972613 # in atomic units
angstrom = 0.1*nm

spcedata = numpy.loadtxt('../../papers/water-saft/figs/sphere-energy-vs-diameter-spce.dat')
dftdata = numpy.loadtxt('../../papers/water-saft/figs/sphere.dat')
diameter_dft = dftdata[:,0]
dftdata[0,1:] = 0 # zero diameter sphere has no 

sphere_area_nm2 = numpy.pi*(diameter_dft*nm + 0.00001)**2

pylab.figure()

pylab.plot(diameter_dft*10, dftdata[:,1]/sphere_area_nm2/mNpermeter, label='$F$')
pylab.plot(diameter_dft*10, (dftdata[:,1] + dftdata[:,4])/sphere_area_nm2/mNpermeter, label='$U$')
pylab.plot(diameter_dft*10, -dftdata[:,4]/sphere_area_nm2/mNpermeter, label='$-TS$')

pylab.axhline(0, color='k', linestyle='-')
pylab.axhline(72, color='k', linestyle=':')
#pylab.xlim(0, 20)
pylab.ylabel('energy per solvent accessible surface area (mN/m)')
pylab.xlabel('sphere diameter ($\AA$)')
pylab.legend(loc='best')

#'figs/sphere.dat' u ($1/2):(-$5/(pi*($1+0.00001)*($1+0.00001)*nm*nm)/mNpermeter) title '-TS' with l ls 3, \
#'figs/sphere.dat' u ($1/2):(($2+$5)/(pi*($1+0.00001)*($1+0.00001)*nm*nm)/mNpermeter) title 'U' with l ls 4, \
#'figs/sphere.dat' u ($1/2):($2/(pi*($1+0.00001)*($1+0.00001)*nm*nm)/mNpermeter) title 'F' with lp ls 1, \

pylab.savefig('figs/sphere-energy.pdf')
