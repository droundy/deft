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
colors = ['k', 'b', 'g', 'r']
plots = ['mc', 'this-work', 'fischer', 'sokolowski'] # , 'gloor'
titles = ['Monte Carlo', 'this work', 'Fischer et al.', 'Sokolowski et al.'] # , 'gloor'
dx = 0.1
sigma = 2.0
contact_delta = 0.1
############################
rpath = sigma + contact_delta

able_to_read_file = True

# Set the max parameters for plotting.
zmax = 9.0
zmin = 1.0
rmax = 4.1
############################

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])
#arg ff = [0.1, 0.2, 0.3, 0.4]

def read_triplet_path(ff, fun):
  if fun == 'mc':
    data = loadtxt("figs/mc/triplet/tripletMC-%03.1f-path2-trimmed.dat" % ff)
    data[:,0]-=4.995
  else:
    # input: "figs/triplet-path-inbetween-*-%4.2f.dat" % (ff)
    filename = "figs/triplet-path-inbetween-%s-%4.2f.dat" % (fun, ff)
    if os.path.isfile(filename) == False:
      data = loadtxt("figs/mc/triplet/tripletMC-%03.1f-path2-trimmed.dat" % ff)
    else:
      data = loadtxt(filename)
      z = data[:,2]
      x = data[:,3]
      zcontact = z.min()
      x[z == zcontact] *= -1
      data[:,3] = x
  return data[:,0:4]

def read_triplet(ff, fun):
  if fun == 'mc':
    return loadtxt("figs/mc/triplet/tripletMC-%3.1f-04.20-trimmed.dat" % ff)

def read_gr(ff):
  return loadtxt("figs/gr-%04.2f.dat" % ff)

fig = figure(figsize=(10,5))

xplot = fig.add_subplot(1,2,2)
zplot = xplot.twiny()
twod_plot = fig.add_subplot(1,2,1)

xlow = -6
xhigh = 8
xplot.set_xlim(xlow, xhigh)
zplot.set_xlim(xlow + rpath, xhigh + rpath)

xticks = [-6, -4, -2, 0]
zticks = [4, 6, 8, 10]

plotticks = xticks + [z-rpath for z in zticks]
xplot.set_xticks(plotticks)
xplot.set_xticklabels(xticks + zticks)
zplot.set_xticks([])


zplot.axvline(x=rpath, color='k')
zplot.axvline(x=3*rpath, color='k')

xloc = .620
zloc = .800
figtext(xloc, .04, r"$\underbrace{\hspace{9.0em}}$", horizontalalignment='center')
figtext(xloc, .01, r"$x$", horizontalalignment='center')
figtext(zloc, .04, r"$\underbrace{\hspace{12.9em}}$", horizontalalignment='center')
figtext(zloc, .01, r"$z$", horizontalalignment='center')

xmin = 1.0
xmax = 9.0
ymax = 6.0

twod_plot.set_xlim(zmin, zmax)
twod_plot.set_ylim(-rmax, rmax)

fig.subplots_adjust(hspace=0.001)

for i in range(len(plots)):
    g3_path = read_triplet_path(ff, plots[i])
    if able_to_read_file == False:
        plot(arange(0,10,1), [0]*10, 'k')
        suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
        savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
        savefig(savedfilename)
        exit(0)
    x = g3_path[:,3]
    z = g3_path[:,2]
    g = g3_path[:,1]
    zcontact = z.min()

    if plots[i] == 'fischer':
      # Fischer et al only predict pair distribution function in
      # contact.  We do this using "&" below which means "and".
      incontact = x**2 + (z-2*rpath)**2 < (rpath + 0.01)**2
      zplot.plot(z[incontact],g[incontact], label=titles[i], color=colors[i])
    else:
      xplot.plot(x[z==zcontact],g[z==zcontact], label=titles[i], color=colors[i])
      zplot.plot(z[z>zcontact],g[z>zcontact], label=titles[i], color=colors[i])


zplot.set_ylabel(r'$g^{(2)}(\left< 0,0,0\right>,\mathbf{r}_2)$')
zplot.legend(loc=1, ncol=2, bbox_to_anchor=(1.2, 1.05), fontsize = 8)

twod_plot.set_aspect('equal')
g2mc = read_triplet(ff, 'mc')
rbins = round(2*rmax/dx) + 2
zminbin = round(zmin/dx)
zmaxbin = round(zmax/dx) + 1
zbins = zmaxbin - zminbin
g22 = zeros((rbins, zbins))
g22[rbins/2:rbins, 0:zbins] = g2mc[:rbins/2, zminbin:zmaxbin]
g22[:rbins/2, 0:zbins] = flipud(g2mc[:rbins/2, zminbin:zmaxbin])
g2mc = g22
gmax = g2mc.max()
dx = 0.1

r = arange(0-rmax, rmax+2*dx, dx)
z = arange(zmin, zmax+dx, dx)
Z, R = meshgrid(z, r)

levels = linspace(0, gmax, gmax*100)
gr = read_gr(ff)
ginf = interp(2*rpath, gr[:,0], gr[:,1]/ff)
xlo = 0.85*ginf/gmax
xhi = 1.15*ginf/gmax
xwhite = 1.0*ginf/gmax
xhier = (1 + xhi)/2.0

cdict = {'red':   [(0.0,  0.0, 0.0),
                   (xlo,  1.0, 1.0),
                   (xwhite,  1.0, 1.0),
                   (xhi,  0.0, 0.0),
                   (xhier,0.0, 0.0),
                   (1.0,  1.0, 1.0)],

         'green': [(0.0, 0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (xwhite, 1.0, 1.0),
                   (xhi, 0.0, 0.0),
                   (xhier,1.0, 1.0),
                   (1.0, 1.0, 1.0)],

         'blue':  [(0.0,  0.0, 0.0),
                   (xlo,  0.1, 0.1),
                   (xwhite,  1.0, 1.0),
                   (xhi,  1.0, 1.0),
                   (xhier,0.0, 0.0),
                   (1.0,  0.0, 0.0)]}
cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)

CS = pcolormesh(Z, R, g2mc, vmax=gmax, vmin=0, cmap=cmap)

sphere0 = Circle((0, 0), 1, color='slategray')
sphere1 = Circle((2*rpath, 0), 1, color='slategray')
twod_plot.add_artist(sphere0)
twod_plot.add_artist(sphere1)

myticks = arange(0, floor(2.0*gmax)/2 + 0.1, 0.5)
colorbar(CS, extend='neither', ticks=myticks)
twod_plot.set_ylabel('$x_2$');
twod_plot.set_xlabel('$z_2$');

# Here we plot the paths on the 2d plot.  The mc plot should align
# with the dft one.
g3_path = read_triplet_path(ff, 'mc')
xmc = g3_path[:,3]
zmc = g3_path[:,2]
g3_path = read_triplet_path(ff, 'this-work')
xdft = g3_path[:,3]
zdft = g3_path[:,2]
plot(zmc,xmc, 'w-', linewidth=3)
plot(zmc,xmc, 'k--', linewidth=3)
plot(zdft,xdft, 'c--', linewidth=3)

Ax = -3.9
Az = rpath
Bx = 0
Bz = rpath
Cx = rpath
Cz = rpath*2
Dx = 0
Dz = rpath*3
Ex = 0
Ez = 8

hw = 4 # headwidth of arrows

g3nice = read_triplet_path(ff, 'this-work')
znice = g3nice[:,2]
xnice = g3nice[:,3]
gnice = g3nice[:,1]
znicecontact = znice.min()
def g3pathfunction_x(x):
    return interp(x, flipud(xnice[znice == znicecontact]), flipud(gnice[znice == znicecontact]))
def g3pathfunction_z(z):
    return interp(z, flipud(znice), flipud(gnice))

# Annotations on 2d plot
annotate('$A$', xy=(Az, Ax), xytext=(1.2,-3.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$B$', xy=(Bz, Bx), xytext=(1.5,1.8), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$C$', xy=(Cz, Cx), xytext=(rpath*2,3.0), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$D$', xy=(Dz, Dx), xytext=(6.7,0.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
annotate('$E$', xy=(Ez, Ex), xytext=(8.5,0.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))


# Annotations on 1d plot
xplot.annotate('$A$', xy=(Ax, g3pathfunction_x(Ax)),
               xytext=(Ax-0.5,0.7),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
xplot.annotate('$B$', xy=(Bx,g3pathfunction_x(Bx)),
               xytext=(Bx+1, g3pathfunction_x(Bx)),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$C$', xy=(Cz,g3pathfunction_z(Cz)),
               xytext=(Cz,g3pathfunction_z(Cz)-0.5),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$D$', xy=(Dz,g3pathfunction_z(Dz)),
               xytext=(Dz+1,g3pathfunction_z(Dz)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('$E$', xy=(Ez,g3pathfunction_z(Ez)),
               xytext=(Ez-0.5,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))


twod_plot.set_title(r'$g^{(3)}(\left< 0,0,0\right>,\left< 0,0,2\sigma\right>,\mathbf{r}_2)$ at $\eta = %g$' % ff)
savefig("figs/triplet-correlation-pretty-inbetween-%d.pdf" % (int(ff*10)))
show()

