#!/usr/bin/python

# We need the following two lines in order for matplotlib to work
# without access to an X server.
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
#import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

#matplotlib.use('Agg')

import pylab, numpy, sys, math

if len(sys.argv) != 3:
    print("Usage:  " + sys.argv[0] + " wb-twodim-filename.dat out-filename.pdf")
    exit(1)

dftdata = numpy.loadtxt(sys.argv[1])

# >>> import matplotlib
# >>> matplotlib.use('Agg')
# >>> import matplotlib.pyplot as plt
# >>> fig=plt.figure()
# >>> fig.save_fig('test.png')

dft_len = len(dftdata[:,0])
dft_dz = dftdata[2,1] - dftdata[1,1]
n0 = dftdata[:,7]
nA = dftdata[:,9]

#step = 0.04
maxval = 1.0
fig = plt.figure()
ax = Axes3D(fig)

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D


# create supporting points in polar coordinates
#r = numpy.linspace(0,1.25,50)
#p = numpy.linspace(0,2*np.pi,50)
#R,P = np.meshgrid(r,p)
# transform them to cartesian system
#X,Y = R*np.cos(P),R*np.sin(P)

dft_len = len(dftdata[:,0])
Nz = int(math.sqrt(dft_len))
#dft_dr = dftdata[2,0] - dftdata[1,0]
z = range(Nz)
y = range(Nz)
yaxis,zaxis = numpy.meshgrid(y,z)
Z,Y = dftdata[Nz*zaxis+yaxis,0],dftdata[Nz*zaxis+yaxis,1]
value = dftdata[Nz*zaxis+yaxis,2]

#Y,Z,value = dftdata[:,0], dftdata[:,1], dftdata[:,2]*4*numpy.pi/3
# n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"r+-",label='wallsDFT $n$')
# n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
# n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")
# pylab.xlim(0,12)
# pylab.ylim(-5,10)

#Z = ((R**2 - 1)**2)
ax.plot_surface(Z, Y, value, rstride=1, cstride=1, cmap=cm.jet)
ax.set_zlim3d(-0.04, 0.04)
ax.set_xlabel(r'$Z position$')
ax.set_ylabel(r'$Y position$')
ax.set_zlabel(r'$Value')
#fig.save_fig('test-particle.pdf')
plt.show()



#pylab.figure(figsize=(8, 8))
#pylab.subplots_adjust(hspace=0.001)

#n_plt = pylab.subplot(1,1,1)
# n_plt.plot(dftdata[:,0],dftdata[:,1]*4*numpy.pi/3,"r+-",label='wallsDFT $n$')
# n_plt.plot(dftdata[:,0],n0*4*numpy.pi/3,"c--",label="$n_0$")
# n_plt.plot(dftdata[:,0],nA*4*numpy.pi/3,"m-.",label="$n_A$")
# pylab.xlim(0,12)
# pylab.ylim(-5,10)
# # if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.45) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.35)):
# #     pylab.ylim(0.1,0.8)
# # if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.35) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.25)):
# #     pylab.ylim(0.17,0.43)
# # if ((mcdata[int(n/2),1]*4*numpy.pi/3<0.15) & (mcdata[int(n/2),1]*4*numpy.pi/3>0.05)):
# #     pylab.ylim(0.04,0.19)
# # pylab.ylim(2.9,3.1)
# pylab.xlabel("position")

# pylab.ylabel("filling fraction")

# pylab.legend(loc='lower left', ncol=2).get_frame().set_alpha(0.5)

# pylab.savefig(sys.argv[2])

# pylab.show()


