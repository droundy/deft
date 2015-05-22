#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.lines as mlines
import os
import numpy

figure()
X = loadtxt('radial-sw-1.00-1.30-0.30-X.dat')
Y = loadtxt('radial-sw-1.00-1.30-0.30-Y.dat')
gradial = loadtxt('radial-sw-1.00-1.30-0.30-n.dat')
contourf(X,Y,gradial,100)
colorbar()
savefig('contour-square-well.pdf')

figure()
data = loadtxt('radial-sw-1.00-1.30-0.30.dat')
plot(data[:,0], data[:,1])
xlabel('r')
ylabel('g(r)')
savefig('radial-square-well.pdf')
show()
