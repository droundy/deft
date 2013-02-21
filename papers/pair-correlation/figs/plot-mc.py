#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
#matplotlib.use('Agg')
import pylab, numpy, sys
import os.path
import matplotlib.colors as mcolors
from matplotlib.widgets import Slider, RadioButtons

ff = 0.3
z0 = 0.05
if (len(sys.argv) > 1):
    ff = float(sys.argv[1])
if (len(sys.argv) > 2):
    z0 = float(sys.argv[2])

zmax = 20
rmax = 10

def read_walls(ff, z0):
    filename = "mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
    print 'Using', filename
    if (os.path.isfile(filename) == False):
        print "File does not exist. Try different values for ff and z0, or leave them blank to use defaults, or generate more monte carlo data."
        sys.exit(1)
    data = numpy.loadtxt(filename)
    return data

g2 = read_walls(ff, z0)

#g2/= 2.552

zbins = len(g2[0,:])
rbins = len(g2[:,0])
dr = rmax/rbins
dz = zmax/zbins

r = numpy.arange(0, rmax, dr)
z = numpy.arange(0, zmax, dz)
Z, R = numpy.meshgrid(z, r)

fig = pylab.figure(1)
ax = pylab.subplot(111)
print g2.max()

levels = 100
CS = pylab.contourf(Z, R, g2, levels)
CS2 = pylab.contourf(Z, -R, g2, levels)

CB = pylab.colorbar(CS)
pylab.axes().set_aspect('equal')

pylab.title('$g^{(2)}(z_0, z_1, r_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
pylab.xlabel("z")
pylab.ylabel("r")

# slider code
z0ax = pylab.axes([0.25, 0, 0.5, 0.05], axisbg='slategray')
z0_slider = Slider(z0ax, '$z_0$', 0.0, 9.90, valinit = 0.0)

def update(val):
  global z0, ax
  ax.collections = []
  z0 = z0_slider.val - z0_slider.val%0.10 + 0.05
  g2 = read_walls(ff, z0)
  CS = ax.contourf(Z, R, g2, levels)
  CS2 = ax.contourf(Z, -R, g2, levels)
  CB.set_clim(vmax = g2.max())
  CB.draw_all()
  ax.set_title('$g^{(2)}(z_0, z_1, r_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
  pylab.draw()

z0_slider.on_changed(update)

# radio button code
ffax = pylab.axes([0.05, 0.3, 0.1, 0.4], axisbg='slategray')
ffbutt = RadioButtons(ffax, (0.1, 0.2, 0.3, 0.4, 0.5), active=2)

def updateff(label):
  global ff, ax, CB
  ax.collections = []
  ff = float(label)
  g2 = read_walls(ff, z0)
  CS = ax.contourf(Z, R, g2, levels)
  CS2 = ax.contourf(Z, -R, g2, levels)
  CB.set_clim(vmax = g2.max())
  CB.draw_all()
  ax.set_title('$g^{(2)}(z_0, z_1, r_1)$, $z_0 = %g$, $ff = %g$' %(z0, ff))
  pylab.draw()

ffbutt.on_clicked(updateff)

pylab.show()
