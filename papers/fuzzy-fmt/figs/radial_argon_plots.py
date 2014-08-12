#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *

figure()
data = loadtxt('figs/YarnellArgon85K.dat')
plot(data[:,0],data[:,1])
savefig('figs/Argon-vapor_pressure-85K.pdf')

figure()
data2 = loadtxt('figs/EggertArgon0.6GPaNN.dat')
plot(data2[:,0],data2[:,1])
savefig('figs/Argon-0_6GPa-RT.pdf')

figure()
data3 = loadtxt('figs/EggertArgon1.1GPaNN.dat')
plot(data3[:,0],data3[:,1])
savefig('figs/Argon-1_1GPa-RT.pdf')
