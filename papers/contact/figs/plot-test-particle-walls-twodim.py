#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt

#contourf plots

import pylab, numpy, sys, math

if len(sys.argv) != 3:
    print("Usage:  " + sys.argv[0] + " wb-twodim-filename.dat out-filename.pdf")
    exit(1)

dftdata = numpy.loadtxt(sys.argv[1])

dft_len = len(dftdata[:,0])
dft_dz = dftdata[2,1] - dftdata[1,1]
n0 = dftdata[:,7]
nA = dftdata[:,9]

#step = 0.04
maxval = 1.0
fig = plt.figure()
#ax = Axes3D(fig)

dft_len = len(dftdata[:,0])
Nz = int(math.sqrt(dft_len))
z = range(Nz)
y = range(Nz)
yaxis,zaxis = numpy.meshgrid(y,z)
Z,Y = dftdata[Nz*zaxis+yaxis,0],dftdata[Nz*zaxis+yaxis,1]
value = dftdata[Nz*zaxis+yaxis,2]

nx = len(dftdata[:,0][dftdata[:,0] == 0])
ny = len(dftdata[:,0])/nx
print 'Nx is', nx
print 'len is', len(dftdata[:,0])
print 'sqr nx is', nx*nx
X = numpy.reshape(dftdata[:,0], (ny,nx))
Y = numpy.reshape(dftdata[:,1], (ny,nx))
n = numpy.reshape(dftdata[:,2], (ny,nx))
pylab.contourf(X, Y, n, 100)
pylab.axes().set_aspect('equal')
pylab.xlabel(r'$Z position$')
pylab.ylabel(r'$Y position$')
plt.show()




