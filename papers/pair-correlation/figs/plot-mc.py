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
#from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot1d():
  global g2, ax
  zvals = len(g2[0][0,:])
  rvals = len(g2[0][:,0])
  if (abs(theta) < numpy.arctan(1/2)):
    z1, r1 = (zmax-zmax/zvals), (zmax-zmax/zvals)*numpy.tan(theta)
  else:
    z1, r1 = (rmax-rmax/rvals)/numpy.tan(theta), (rmax-rmax/rvals)

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

colors = ['r', 'k', 'b', 'cyan']
plots = ['mc', 'simple', 'nA', 'py']
numplots = 4

def plot():
  global ax, CS
  i = 0
  while i < numplots:
    ax[i].collections = []
    g2[i] = read_walls(ff, z0, plots[i])
    CS = ax[i].contourf(Z, R, g2[i], levels, cmap=cmap, extend='both')
    CS2 = ax[i].contourf(Z, -R, g2[i], levels, cmap=cmap, extend='both')
    CS.cmap.set_over('k')
    CS2.cmap.set_over('k')
    ax[i].set_title('%s, $z_0 = %g$, $ff = %g$' %(plots[i], z0, ff))
    i += 1
  plot1d()
  pylab.draw()

ff = 0.3
z0 = 0.05

zmax = 20
rmax = 10

def read_walls(ff, z0, fun):
  if fun == 'mc':
    filename = "mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
  else:
    filename = "walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
  print 'Using', filename
  if (os.path.isfile(filename) == False):
    print "File does not exist."
    sys.exit(1)
  data = numpy.loadtxt(filename)
  return data

# COLORS
cdict = {'red':  ((0.0, 0.0, 0.0),
                  (0.25,0.56, 0.2),
                  (0.5, 1.0, 1.0),
                  (0.75,1.0, 1.0),
                  (1.0, 1.0, 0.0)),

        'green': ((0.0, 0.0, 0.0),
                  (0.25,1.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (0.75,0.0, 0.0),
                  (1.0, 1.0, 0.0)),

        'blue':  ((0.0, 0.0, 0.0),
                  (0.25,0.6, 1.0),
                  (0.5, 1.0, 1.0),
                  (0.75,0.0, 0.0),
                  (1.0, 1.0, 0.0))
        }
map = LinearSegmentedColormap('map', cdict)
pylab.register_cmap(cmap = map)
cmap = pylab.get_cmap('map')

levels = numpy.linspace(0, 4, 49)

g2 = [0]*numplots
g2[0] = read_walls(ff, z0, 'mc')
theta = 0
zbins = len(g2[0][0,:])
rbins = len(g2[0][:,0])
dr = rmax/rbins
dz = zmax/zbins

r = numpy.arange(0, rmax, dr)
z = numpy.arange(0, zmax, dz)
Z, R = numpy.meshgrid(z, r)

fig = pylab.figure(1)
ax=[0]*5
left = .05
right = .95
bottom = .1
top = .96
hspace = .02
vspace = .04
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
z0ax = pylab.axes([0.25, 0.01, 0.5, 0.025], axisbg='slategray')
z0_slider = Slider(z0ax, 'z$_0$', 0.0, 9.90, valinit = 0.0)

def update(val):
  global z0
  z0 = z0_slider.val - z0_slider.val%0.10 + 0.05
  plot()
z0_slider.on_changed(update)

# angle slider
angax = pylab.axes([0.25, 0.035, 0.5, 0.025], axisbg='slategray')
angslider = Slider(angax, 'theta', 0, numpy.pi/2, valinit = numpy.pi/2)

def upangle(val):
  global theta
  theta = numpy.pi/2 - angslider.val
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

pylab.show()
