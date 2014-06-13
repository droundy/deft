#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
#matplotlib.use('Agg')
import pylab, numpy, sys, scipy.ndimage
import os.path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider, RadioButtons
#from mpl_toolkits.axes_grid1 import make_axes_locatab

at_wall = True
# these are the things to set
colors = ['k', 'b', 'r', 'g']
plots = ['mc', 'this-work', 'fischer', 'gloor']
numplots = 4
if at_wall == True:
  colors = ['k', 'b', 'r', 'g']
  #plots = ['mc', 'this-work', 'fischer', 'sphere-dft']
  plots = ['mc', 'this-work', 'this-work-mc', 'sokolowski']
  numplots = 4

dx = 0.1
############################
z0 = 0.05
theta = 0
ff = 0.3

def plot1d():
  global g2, ax
  zvals = len(g2[0][0,:])
  rvals = len(g2[0][:,0])
  if theta < numpy.arctan(1/2):
    z1 = zmax
    r1 = (z1-z0)*numpy.tan(theta)
  elif theta < numpy.pi - numpy.arctan(1/2):
    r1 = rmax
    z1 = r1/numpy.tan(theta) + z0
  else:
    z1 = -zmax
    r1 = (z1-z0)*numpy.tan(theta)
  rlen = numpy.sqrt((z1-z0)**2 + r1**2)
  z0coord = zvals*(z0-0.05)/zmax
  z1coord = zvals*(z1)/zmax
  r0coord = 0
  r1coord = rvals*(r1)//rmax
  y, x = numpy.linspace(z0coord, z1coord, num), numpy.linspace(r0coord, r1coord, num)

  i = 0
  while i < numplots:
    zi = scipy.ndimage.map_coordinates(g2[i], numpy.vstack((x,y)), order = 1)
    angline[i].set_data([z0, z1], [0, r1])
    gslice[i].set_data(numpy.linspace(0, rlen, num), zi)
    i += 1

def plot():
  global ax, CS
  i = 0
  while i < numplots:
    ax[i].collections = []
    g2[i] = read_walls(ff, z0, plots[i])

    rmax = len(g2[i][:,0])*dx
    zmax = len(g2[i][0,:])*dx
    print "zmax len is = %d and zmax is = %g" % (len(g2[i][0,:]), zmax)

    r = numpy.arange(0, rmax, dx)
    z = numpy.arange(0, zmax, dx)
    Z, R = numpy.meshgrid(z, r)

    CS = ax[i].contourf(Z, R, g2[i], levels, cmap=cmap)
    CS2 = ax[i].contourf(Z, -R, g2[i], levels, cmap=cmap)
    #CS.cmap.set_over('k')
    #CS2.cmap.set_over('k')
    ax[i].set_title('%s, $z_0 = %g$, $ff = %g$' %(plots[i], z0, ff))
    i += 1
  plot1d()
  pylab.draw()

def read_walls(ff, z0, fun):
  if fun == 'mc':
    filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
  elif fun == 'sphere-dft':
    filename = "figs/wallsWB-with-sphere-%04.2f.dat" % ff
  else:
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
  print 'Using', filename
  if (os.path.isfile(filename) == False):
    print "File does not exist:", filename
    sys.exit(1)
  data = numpy.loadtxt(filename)
  return data

gmax = 4
xlo = 0.85/gmax
xhi = 1.15/gmax
xhier = (1+xhi)/2.0

cdict = {'red':   [(0.0,  0.0, 0.0),
                   (xlo,  1.0, 1.0),
                   (1.0/gmax,  1.0, 1.0),
                   (xhi,  0.0, 0.0),
                   (xhier,0.0, 0.0),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0, 0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (1.0/gmax, 1.0, 1.0),
                   (xhi, 0.0, 0.0),
                   (xhier,1.0, 1.0),
                   (1.0, 1.0, 1.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (1.0/gmax,  1.0, 1.0),
                   (xhi,  1.0, 1.0),
                   (xhier,0.0, 0.0),
                   (1.0,  0.0, 0.0)]}
cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)
levels = numpy.linspace(0, gmax, gmax*100)

g2 = [0]*numplots
g2[0] = read_walls(ff, z0, 'mc')

rmax = len(g2[0][:,0])*dx
zmax = len(g2[0][0,:])*dx

r = numpy.arange(0, rmax, dx)
z = numpy.arange(0, zmax, dx)
Z, R = numpy.meshgrid(z, r)

fig = pylab.figure(1)
ax=[0]*(numplots+1)
left = .05
right = .95
bottom = .1
top = .96
hspace = .02
vspace = .1
width = (right - left - 2*hspace)/3
height = (top - bottom - vspace)/2

ax[0] = pylab.axes([left, top-height, width, height])
ax[1] = pylab.axes([left, bottom, width, height])
ax[2] = pylab.axes([left+width+hspace, bottom, width, height])
ax[3] = pylab.axes([left+width*2+hspace*2, bottom, width, height])
ax[4] = pylab.axes([left+hspace+width, top-height, width*2+hspace, height])

i = 0
while i < numplots:
  ax[i].set_aspect('equal')
  ax[i].set_xlim(0, 20)
  ax[i].set_ylim(-10,10)
  i += 1

num = 1000

angline = [0]*numplots
gslice = [0]*numplots
i = 0
while i < numplots:
  angline[i], = ax[i].plot([0, zmax], [0, rmax], 'k')
  gslice[i], = ax[4].plot(numpy.zeros(num), numpy.zeros(num), colors[i])
  i += 1

ax[4].set_xlim(2, 6.5)
ax[4].set_ylim(0, 4)

pylab.legend(plots)
ax[4].axhline(y=1, linestyle='--', color='slategray')

plot()
#pylab.tight_layout()

# colorbar
#divider = make_axes_locatable(ax[0])
#cax = divider.append_axes("right", "5%", pad="3%")
#ticks = numpy.linspace(0, 5, 11)
#CB = pylab.colorbar(CS, cax=cax)
#CB.set_ticks(ticks)

# slider
z0_valinit = 5.0
z0ax = pylab.axes([0.25, 0.01, 0.5, 0.025], axisbg='slategray')
z0_slider = Slider(z0ax, 'z$_0$', 0.0, 9.90, valinit = z0_valinit)

def update(val):
  global z0
  z0 = z0_slider.val - z0_slider.val%0.10 + 0.05
  plot()
z0_slider.on_changed(update)

# angle slider
angle_valinit = numpy.pi/3
angax = pylab.axes([0.25, 0.035, 0.5, 0.025], axisbg='slategray')
angslider = Slider(angax, 'theta', 0, numpy.pi, valinit = 0)


def upangle(val):
  global theta
  theta = angslider.val
  plot1d()

angslider.on_changed(upangle)

# radio buttons
ffax = pylab.axes([0.01, 0.1, 0.04, 0.2], axisbg='slategray')
ffbutt = RadioButtons(ffax, ('.1', '.2', '.3', '.4', '.5'), active=2)

def updateff(label):
  global ff
  ff = float(label)
  plot()
ffbutt.on_clicked(updateff)

savedfilename = "pair-correlation-" + str(int(ff*10)) + "-" + str(int(z0_valinit*10)) + "-" + str(int(10*angle_valinit/numpy.pi)) + ".pdf"
#pylab.savefig(savedfilename)
pylab.show()
