#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *

sigma = 0.3506 #nm

figure()
data = loadtxt('figs/YarnellArgon85K.dat')
n =0.02125 # Angstrom^-3
nsig_3 = n*(sigma*10)**3
plot(data[:,0],data[:,1]*nsig_3)
savefig('figs/Argon-vapor_pressure-85K.pdf')

figure()
data2 = loadtxt('figs/EggertArgon0.6GPaNN.dat')
n = 24.23 #atoms/nm^3
nsig_3 = n*sigma**3
plot(data2[:,0],data2[:,1]*nsig_3)
savefig('figs/Argon-0_6GPa-RT.pdf')

figure()
data3 = loadtxt('figs/EggertArgon1.1GPaNN.dat')
n = 27.74 #atoms/nm^3
nsig_3 = n*sigma**3
plot(data3[:,0],(data3[:,1]-1)*nsig_3)
savefig('figs/Argon-1_1GPa-RT.pdf')
