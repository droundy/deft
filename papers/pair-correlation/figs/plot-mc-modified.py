#!/usr/bin/python
# Produces a contour plot for g^(2) (z0, z1, x1).
# Takes z0 and filling fraction as arguments.

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import pylab, numpy, sys, scipy.ndimage
import os.path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider, RadioButtons
#from mpl_toolkits.axes_grid1 import make_axes_locatable

colors = ['k', 'b', 'r', 'g']
plots = ['mc', 'this-work', 'fischer', 'gloor']
numplots = 4
dx = 0.1
able_to_read_file = True
z0 = 0.95
theta = numpy.pi/2
ff = 0.3

if len(sys.argv) != 4:
    print("Usage:  " + sys.argv[0] + " ff z0 theta/pi")
    exit(1)
ff = float(sys.argv[1])
z0 = float(sys.argv[2])
#theta_div_by_pi = float(sys.argv[2])
theta = numpy.pi*float(sys.argv[3])
print theta/numpy.pi

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
    ax[4].set_title('$z_0 = %g$, $ff = %g$, $\\theta = %g$' %(z0, ff, theta))

    i = 0
    ax[4].set_title('%s, $z_0 = %g$, $ff = %g$ $theta = %g\pi$' %(plots[i], z0, ff, theta/numpy.pi))
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
        global able_to_read_file
        if able_to_read_file == False:
            break

        rmax = len(g2[i][:,0])*dx
        zmax = len(g2[i][0,:])*dx

        r = numpy.arange(0, rmax, dx)
        z = numpy.arange(0, zmax, dx)
        Z, R = numpy.meshgrid(z, r)

        CS = ax[i].contourf(Z, R, g2[i], levels, cmap=cmap, extend='both')
        CS2 = ax[i].contourf(Z, -R, g2[i], levels, cmap=cmap, extend='both')
        CS.cmap.set_over('k')
        CS2.cmap.set_over('k')
        ax[i].set_title('%s, $z_0 = %g$, $ff = %g$' %(plots[i], z0, ff))
        i += 1
    plot1d()
    pylab.draw()

def read_walls(ff, z0, fun):
    if fun == 'mc':
        filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            global able_to_read_file
            able_to_read_file = False
            return 0
    else:
        filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            global able_to_read_file
            able_to_read_file = False
            return 0
    print 'Using', filename
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

if able_to_read_file == False:
    matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
    matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
    savedfilename = "figs/pair-correlation-" + str(int(ff*10)) + "-" + str(int(z0*100)) + "-" + str(int(10*theta/numpy.pi)) + ".pdf"
    pylab.savefig(savedfilename)
    exit(0)

rmax = len(g2[0][:,0])*dx
zmax = len(g2[0][0,:])*dx

print rmax, zmax


zbins = len(g2[0][0,:])
rbins = len(g2[0][:,0])
dr = rmax/rbins
dz = zmax/zbins

r = numpy.arange(0, rmax, dx)
z = numpy.arange(0, zmax, dx)
Z, R = numpy.meshgrid(z, r)

#fig = pylab.figure(1)
ax=[0]*5
left = .05
right = .95
bottom = .1
top = .96
hspace = .05
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
if able_to_read_file == False:
    matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
    matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
    savedfilename = "figs/pair-correlation-" + str(int(ff*10)) + "-" + str(int(z0*100)) + "-" + str(int(10*theta/numpy.pi)) + ".pdf"
    pylab.savefig(savedfilename)
    exit(0)

plot()


savedfilename = "figs/pair-correlation-" + str(int(ff*10)) + "-" + str(int(z0*100)) + "-" + str(int(10*theta/numpy.pi)) + ".pdf"
pylab.savefig(savedfilename)
pylab.show()


#pylab.tight_layout()

# colorbar
#divider = make_axes_locatable(ax[0])
#cax = divider.append_axes("right", "5%", pad="3%")
#ticks = numpy.linspace(0, 5, 11)
#CB = pylab.colorbar(CS, cax=cax)
#CB.set_ticks(ticks)

# slider
# z0_valinit = 5.0
# z0ax = pylab.axes([0.25, 0.01, 0.5, 0.025], axisbg='slategray')
# z0_slider = Slider(z0ax, 'z$_0$', 0.0, 9.90, valinit = z0_valinit)

# def update(val):
#   global z0
#   z0 = z0_slider.val - z0_slider.val%0.10 + 0.05
#   plot()
# z0_slider.on_changed(update)

# # angle slider
# angle_valinit = numpy.pi/3
# angax = pylab.axes([0.25, 0.035, 0.5, 0.025], axisbg='slategray')
# angslider = Slider(angax, 'theta', 0, numpy.pi/2, valinit = angle_valinit)


# def upangle(val):
#   global theta
#   theta = numpy.pi/2 - angslider.val
#   plot1d()

# angslider.on_changed(upangle)

# radio buttons
# ffax = pylab.axes([0.01, 0.1, 0.04, 0.2], axisbg='slategray')
# ffbutt = RadioButtons(ffax, ('.1', '.2', '.3', '.4', '.5'), active=2)

# def updateff(label):
#   global ff
#   ff = float(label)
#   plot()
# ffbutt.on_clicked(updateff)
