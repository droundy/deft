#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt

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
ax = Axes3D(fig)

dft_len = len(dftdata[:,0])
Nz = int(math.sqrt(dft_len))
z = range(Nz)
y = range(Nz)
yaxis,zaxis = numpy.meshgrid(y,z)
Z,Y = dftdata[Nz*zaxis+yaxis,0],dftdata[Nz*zaxis+yaxis,1]
value = dftdata[Nz*zaxis+yaxis,2]

ax.plot_wireframe(Z, Y, value, rstride=1, cstride=1, cmap=cm.jet)
ax.set_zlim3d(-0.04, 0.04)
ax.set_xlabel(r'$Z position$')
ax.set_ylabel(r'$Y position$')
ax.set_zlabel(r'$Value')
plt.show()




