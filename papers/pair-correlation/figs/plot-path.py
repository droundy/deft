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

# these are the things to set
colors = ['k', 'b', 'g', 'r', 'm']
plots = ['mc', 'this-work', 'sphere-dft', 'fischer', 'sokolowski'] # , 'gloor'
titles = ['Monte Carlo', 'this work', 'test particle', 'Fischer et al.', 'Sokolowski'] # , 'gloor'

dx = 0.1
############################
rpath = 2.005

able_to_read_file = True

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

def read_walls_path(ff, fun):
  if fun == 'mc':
    # input:  "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % (ff)
    filename = "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % ff
  elif fun == 'sphere-dft':
    filename = "figs/wallsWB-with-sphere-path-%1.2f.dat" % 0.3 # ff FIXME others don't exist in repo yet
  else:
    # input: "figs/walls.dat" % ()
    # input: "figs/walls/wallsWB-path-*-pair-%1.2f-0.005.dat" %(ff)
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-0.005.dat" %(fun, ff)
  data = loadtxt(filename)
  if fun == 'mc':
    data[:,0]-=4.995
  else:
    data[:,2]-=3
  return data[:,0:4]

def read_walls_mc(ff):
  # The 0.05 below is deceptive (ask Paho).
  return loadtxt("figs/mc/wallsMC-pair-%1.1f-0.05-trimmed.dat" % ff)

def read_walls_dft(ff, fun):
  if fun == 'sphere-dft':
    filename = "figs/wallsWB-with-sphere-%1.2f-trimmed.dat" % 0.3 # ff FIXME others don't exist in repo yet
  else:
    # input: "figs/walls/wallsWB-*-pair-%1.2f-0.005.dat" %(ff)
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-0.005.dat" %(fun, ff)
  return loadtxt(filename)

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
    g2_path = read_walls_path(ff, plots[i])
    if able_to_read_file == False:
        plot(arange(0,10,1), [0]*10, 'k')
        suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
        savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
        savefig(savedfilename)
        exit(0)
    x = g2_path[:,3]
    z = g2_path[:,2]
    g = g2_path[:,1]
    zcontact = z.min()
    xcontact = 2.0051

    if plots[i] == 'fischer':
      # Fischer et al only predict pair distribution function in
      # contact.  We do this using "&" below which means "and".
      incontact = (x<xcontact) & (z<2)
      zplot.plot(z[incontact],g[incontact], label=titles[i], color=colors[i])
    else:
      xplot.plot(x[z==zcontact],g[z==zcontact], label=titles[i], color=colors[i])
      zplot.plot(z[x<xcontact],g[x<xcontact], label=titles[i], color=colors[i])

zplot.axvline(x=2, color='k')
zplot.axvline(x=0, color='k')


zplot.set_ylabel(r'$g^{(2)}(\left< 0,0,0\right>,\mathbf{r}_2)$')
zplot.legend(loc=1, ncol=2, bbox_to_anchor=(1.25, 1.05), fontsize = 8)

twod_plot.set_aspect('equal')
g2mc = read_walls_mc(ff)
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

Ax = 3.9
Az = 0
Bx = rpath
Bz = 0
Cx = rpath/sqrt(2.0)
Cz = rpath/sqrt(2.0)
Dx = 0
Dz = rpath
Ex = 0
Ez = 4.0

hw = 4 # headwidth of arrows

g2nice = read_walls_path(ff, 'this-work')
def g2pathfunction_x(x):
    return interp(x, flipud(g2nice[:,3]), flipud(g2nice[:,1]))
def g2pathfunction_z(z):
    return interp(z, g2nice[:,2], g2nice[:,1])

# Annotations on 2d plot
annotate('$A$', xy=(Az,Ax), xytext=(1,3), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$B$', xy=(Bz,Bx), xytext=(1,2.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$C$', xy=(Cz,Cx), xytext=(2.3,2.0), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$D$', xy=(Dz,Dx), xytext=(3,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$E$', xy=(Ez,Ex), xytext=(5,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))

# Annotations on 1d plot
xplot.annotate('$A$', xy=(Ax, g2pathfunction_x(Ax)),
               xytext=(Ax+1,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
xplot.annotate('$B$', xy=(Bx,g2pathfunction_x(Bx)),
               xytext=(Bx+1, g2pathfunction_x(Bx)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$C$', xy=(Cz,g2pathfunction_z(Cz)),
               xytext=(Cz,g2pathfunction_z(Cz)-0.5),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$D$', xy=(Dz,g2pathfunction_z(Dz)),
               xytext=(Dz+1,g2pathfunction_z(Dz)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$E$', xy=(Ez,g2pathfunction_z(Ez)),
               xytext=(Ez+0.7,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))


twod_plot.set_title(r'$g^{(2)}(\left< 0,0,0\right>, \left<x_2, 0, z_2\right>)$ at $\eta = %g$' % ff)
savefig("figs/pair-correlation-pretty-%d.pdf" % (int(ff*10)))
show()

