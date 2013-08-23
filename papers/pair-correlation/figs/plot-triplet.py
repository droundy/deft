#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib, scipy.ndimage, os.path
#matplotlib.use('Agg')
from pylab import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider, RadioButtons
from mpl_toolkits.axes_grid1 import make_axes_locatable

# these are the things to set
colors = ['k', 'b', 'r', 'g']
plots = ['this-work', 'fischer', 'gloor']
numplots = 1
dx = 0.1
############################
#z0 = 6.05
z0 = float(sys.argv[1])
theta = 0
ff = 0.3

def plot1d():
  global g2, ax
  zvals = len(g2[0][0,:])
  rvals = len(g2[0][:,0])
  if theta < arctan(1/2):
    z1 = zmax
    r1 = (z1-z0)*tan(theta)
  elif theta < pi - arctan(1/2):
    r1 = rmax
    z1 = r1/tan(theta) + z0
  else:
    z1 = -zmax
    r1 = (z1-z0)*tan(theta)
  rlen = sqrt((z1-z0)**2 + r1**2)
  z0coord = zvals*(z0-0.05)/zmax
  z1coord = zvals*(z1)/zmax
  r0coord = 0
  r1coord = rvals*(r1)//rmax
  y, x = linspace(z0coord, z1coord, num), linspace(r0coord, r1coord, num)

  i = 0
  while i < numplots:
    zi = scipy.ndimage.map_coordinates(g2[i], vstack((x,y)), order = 1)
    angline[i].set_data([z0, z1], [0, r1])
    gslice[i].set_data(linspace(0, rlen, num), zi)
    i += 1

def make_plots():
  global ax, CS
  i = 0
  while i < numplots:
    ax[i].collections = []
    g2[i] = read_walls(ff, z0, plots[i])

    rmax = len(g2[i][:,0])*dx-dx
    zmax = len(g2[i][0,:])*dx-dx

    r = arange(0, rmax, dx)
    z = arange(0, zmax, dx)
    Z, R = meshgrid(z, r)
    print Z.shape
    print R.shape
    print g2[i].shape
    CS = ax[i].contourf(Z, R, g2[i], levels, cmap=cmap)
    CS2 = ax[i].contourf(Z, -R, g2[i], levels, cmap=cmap)
    ax[i].set_title('%s, $z_0 = %g$, $ff = %g$' %(plots[i], z0, ff))
    i += 1
  #plot1d()
  draw()

def read_walls(ff, z0, fun):
  if fun == 'mc':
    filename = "figs/mc/tripletWB-this-w-mc-%03.1f-%05.2f.dat" % (ff, z0)
  else:
    filename = "figs/tripletWB-%s-%1.2f-%1.1f0.dat" %(fun, ff, z0)
  print 'Using', filename
  if (os.path.isfile(filename) == False):
    print "File does not exist:", filename
    return zeros((150,150))
  data = loadtxt(filename)
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
levels = linspace(0, gmax, gmax*20)
levels = linspace(0,4, 100);


g2 = [0]*numplots
g2[0] = read_walls(ff, z0, 'mc')

rmax = len(g2[0][:,0])*dx
zmax = len(g2[0][0,:])*dx

# r = arange(0, rmax, dx)
# z = arange(-zmax, zmax/2, dx)
# Z, R = meshgrid(z, r)

fig = figure(1)
ax=[0]*5
left = .05
right = .95
bottom = .1
top = .96
hspace = .02
vspace = .1
width = (right - left - 2*hspace)/3
height = (top - bottom - vspace)/2

ax[0] = axes([left, top-height, width, height])
ax[1] = axes([left, bottom, width, height])
ax[2] = axes([left+width+hspace, bottom, width, height])
ax[3] = axes([left+width*2+hspace*2, bottom, width, height])
ax[4] = axes([left+hspace+width, top-height, width*2+hspace, height])

i = 0
while i < numplots:
  ax[i].set_aspect('equal')
  ax[i].set_xlim(0, 15)
  ax[i].set_ylim(-7.5,7.5)
  i += 1

num = 1000

angline = [0]*numplots
gslice = [0]*numplots
i = 0
while i < numplots:
  angline[i], = ax[i].plot([0, zmax], [0, rmax], 'k')
  gslice[i], = ax[4].plot(zeros(num), zeros(num), colors[i])
  i += 1

ax[4].set_xlim(2, 6.5)
ax[4].set_ylim(0, 4)

legend(plots)
ax[4].axhline(y=1, linestyle='--', color='slategray')

make_plots()
#tight_layout()

# colorbar
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", "5%", pad="3%")
ticks = linspace(0, 5, 11)
CB = colorbar(CS, cax=cax)
CB.set_ticks(ticks)

# slider
z0_valinit = 2.05
z0ax = axes([0.25, 0.01, 0.5, 0.025], axisbg='slategray')
z0_slider = Slider(z0ax, 'z$_0$', 2.05, 9.95, valinit = z0_valinit)

def update(val):
  global z0
  z0 = z0_slider.val - z0_slider.val%0.10 + 0.05
  make_plots()
z0_slider.on_changed(update)

# angle slider
angle_valinit = 0
angax = axes([0.25, 0.035, 0.5, 0.025], axisbg='slategray')
angslider = Slider(angax, 'theta', 0, pi, valinit = 0)


def upangle(val):
  global theta
  theta = angslider.val
  plot1d()

angslider.on_changed(upangle)

# radio buttons
ffax = axes([0.01, 0.1, 0.04, 0.2], axisbg='slategray')
ffbutt = RadioButtons(ffax, ('.1', '.2', '.3', '.4', '.5'), active=2)

def updateff(label):
  global ff
  ff = float(label)
  make_plots()
ffbutt.on_clicked(updateff)

show()
