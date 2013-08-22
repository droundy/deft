#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if len(sys.argv) < 3 or sys.argv[2] != "show":
  matplotlib.use('Agg')
from pylab import *
import scipy.ndimage
import os.path
import math

from matplotlib import rc
rc('text', usetex=True)

from matplotlib.colors import NoNorm

at_wall = True

# these are the things to set
colors = ['k', 'b', 'g', 'r']
plots = ['mc', 'this-work', 'fischer', 'sokolowski'] # , 'gloor'
titles = ['Monte Carlo', 'this work', 'Fischer et al', 'Sokolowski'] # , 'gloor'
if at_wall:
  plots = ['mc', 'this-work', 'sphere-dft', 'fischer', 'sokolowski'] # , 'gloor'
  titles = ['Monte Carlo', 'this work', 'sphere-dft', 'Fischer et al','Sokolowski'] # , 'gloor'
  colors = ['k', 'b', 'g', 'r','m']
dx = 0.1
############################

able_to_read_file = True
z0 = 0.05

# Set the max parameters for plotting.
zmax = 6
zmin = -.5
rmax = 4.1
############################

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])
#arg ff = [0.1, 0.2, 0.3, 0.4]

def read_walls_path(ff,z0,fun):
  if fun == 'mc':
    # input:  "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % (ff)
    filename = "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % ff
  elif fun == 'sphere-dft':
    filename = "figs/walls/wallsWB-sphere-dft-path-%1.2f.dat" % ff
  else:
    # input: "figs/walls.dat" % ()
    # input: "figs/walls/wallsWB-path-*-pair-%1.2f-*.dat" %(ff)
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
  if (os.path.isfile(filename) == False):
    # Just use walls data if we do not have the MC (need to be careful!)
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-%1.2f.dat" %('this-work', ff, z0)
  data = loadtxt(filename)
  if fun == 'mc':
    data[:,0]-=4.995
  else:
    data[:,2]-=3
  return data[:,0:4]

def read_walls(ff, z0, fun):
  if fun == 'mc':
    # input: "figs/mc/wallsMC-pair-%1.1f-*.dat" % (ff)
    filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
  elif fun == 'sphere-dft':
    filename = "figs/walls/wallsWB-sphere-dft-%1.2f.dat" % ff
  else:
    # input: "figs/walls/wallsWB-*-pair-%1.2f-*.dat" %(ff)
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
  #print 'Using', filename
  if (os.path.isfile(filename) == False):
    # Just use walls data if we do not have the MC (need to be careful!)
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %('this-work', ff, z0)
    title("Using this work instead of MC!")
  data = loadtxt(filename)
  return data

ymax = 4

fig = figure(figsize=(10,5))

xplot = fig.add_subplot(1,2,2)
zplot = xplot.twiny()
#zplot = fig.add_subplot(1,3,3, sharey=xplot)
twod_plot = fig.add_subplot(1,2,1)

xplot.set_xlim(6, -4)
xplot.set_xticks([6, 4, 2, 0, -2, -4])
xplot.set_xticklabels([6, 4, "2 0", 2, 4, 6])
zplot.set_xlim(-4, 6)
zplot.set_xticks([])
#xplot.set_ylim(0)

figtext(.613, .04, r"$\underbrace{\hspace{9em}}$", horizontalalignment='center')
figtext(.613, .01, r"$x$", horizontalalignment='center')
figtext(.796, .04, r"$\underbrace{\hspace{13.1em}}$", horizontalalignment='center')
figtext(.796, .01, r"$z$", horizontalalignment='center')

twod_plot.set_xlim(-0.5, 1.5*ymax)
twod_plot.set_ylim(-ymax, ymax)

fig.subplots_adjust(hspace=0.001)

for i in range(len(plots)):
    g2_path = read_walls_path(ff, z0, plots[i])
    if able_to_read_file == False:
        plot(arange(0,10,1), [0]*10, 'k')
        suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
        savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
        savefig(savedfilename)
        exit(0)
    # FIXME:  once we have sufficient data, we will want to remove this smoothing...
    old = 1.0*g2_path[:,1]
    averaged = 1.0*old
    # averaged[1:] += old[:len(averaged)-1]
    # averaged[:len(averaged)-1] += old[1:]
    # averaged /= 3;

    # Find x-portion:
    z_beg = g2_path[0,2]
    coord = 0
    z_final = z_beg
    print g2_path[:,2]
    print plots[i]
    while z_final == z_beg:
      coord += 1
      z_final = g2_path[coord,2]

    xplot.plot(g2_path[:coord,3],averaged[:coord], label=titles[i], color=colors[i])
    zplot.plot(g2_path[coord-1:,2],averaged[coord-1:], label=titles[i], color=colors[i])

g2nice = read_walls_path(ff, z0, 'this-work')

def g2pathfunction_x(x):
    return interp(x, flipud(g2nice[:,3]), flipud(g2nice[:,1]))
def g2pathfunction_z(z):
    return interp(z, g2nice[:,2], g2nice[:,1])

rA = 3.9
rE = 4.0
rpath = 2.005
xAoff = 3.8
xBoff = 2.0
zCoff = 1.0
zDoff = 2.0
zEoff = 3.8

zplot.axvline(x=2, color='k')
zplot.axvline(x=0, color='k')

hw = 4 # headwidth of arrows

xplot.annotate('$A$', xy=(xAoff, g2pathfunction_x(xAoff)),
               xytext=(xAoff+1,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
xplot.annotate('$B$', xy=(xBoff,g2pathfunction_x(xBoff)),
               xytext=(xBoff+1, g2pathfunction_x(xBoff)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$C$', xy=(zCoff,g2pathfunction_z(zCoff)),
               xytext=(zCoff,g2pathfunction_z(zCoff)-0.5),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$D$', xy=(zDoff,g2pathfunction_z(zDoff)),
               xytext=(zDoff+1,g2pathfunction_z(zDoff)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$E$', xy=(zEoff,g2pathfunction_z(zEoff)),
               xytext=(zEoff+0.7,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))


zplot.set_ylabel(r'$g^{(2)}(\left< 0,0,0\right>,\mathbf{r}_2)$')
zplot.legend(loc=3, ncol=2)

twod_plot.set_aspect('equal')
g2mc = read_walls(ff, z0, 'mc')
rbins = round(2*rmax/dx)
zposbins = round(zmax/dx)
znegbins = round(-zmin/dx)
zbins = zposbins + znegbins
g22 = zeros((rbins, zbins))
g22[rbins/2:rbins, znegbins:zbins] = g2mc[:rbins/2,:zposbins]
g22[:rbins/2, znegbins:zbins] = flipud(g2mc[:rbins/2,:zposbins])
g2mc = g22
gmax = g2mc.max()
dx = 0.1

r = arange(0-rmax, rmax, dx)
z = arange(zmin, zmax, dx)
Z, R = meshgrid(z, r)

levels = linspace(0, gmax, gmax*100)
xlo = 0.85/gmax
xhi = 1.15/gmax
xhier = (1 + xhi)/2.0

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

CS = pcolormesh(Z, R, g2mc, vmax=gmax, vmin=0, cmap=cmap)

myticks = arange(0, floor(2.0*gmax)/2 + 0.1, 0.5)
colorbar(CS, extend='neither', ticks=myticks)
twod_plot.set_ylabel('$x_2$');
twod_plot.set_xlabel('$z_2$');

xs = [0, 0]
ys = [ymax, rpath]
dtheta = pi/80
for theta in arange(0, pi/2 + dtheta/2, dtheta):
    xs.append(rpath*sin(theta))
    ys.append(rpath*cos(theta))
xs.append(2*ymax)

ys.append(0)
plot(xs, ys, 'w-', linewidth=2)
plot(xs, ys, 'k--', linewidth=2)

annotate('$A$', xy=(0,rA), xytext=(1,3), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$B$', xy=(0,rpath), xytext=(1,2.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$C$', xy=(rpath/sqrt(2.0),rpath/sqrt(2.0)), xytext=(2.3,2.0), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$D$', xy=(rpath,0), xytext=(3,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$E$', xy=(rE,0), xytext=(5,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))


twod_plot.set_title(r'$g^{(2)}(\left< 0,0,0\right>, \left<x_2, 0, z_2\right>)$ at $\eta = %g$' % ff)
savefig("figs/pair-correlation-pretty-alt-%d.pdf" % (int(ff*10)))
show()

