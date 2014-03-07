#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if len(sys.argv) < 3 or sys.argv[2] != "show":
  matplotlib.use('Agg')
from pylab import *
import scipy.ndimage
import os.path
import math
import bracket # our handy bracket function
import styles # our preferred line styles
import matplotlib.patheffects
from matplotlib import rc
rc('text', usetex=True)

from matplotlib.colors import NoNorm

# these are the things to set
plots = ['mc', 'this-work', 'this-work-mc', 'sokolowski', 'fischer'] # , 'gloor'
dx = 0.1
sigma = 2.0
contact_delta = 0.1
############################
rpath = sigma + contact_delta
center = rpath

able_to_read_file = True

# Set the max parameters for plotting.
zmax = 9.0
zmin = -4.0
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
  return data[:,0:4]

def read_triplet_back(ff, fun):
  if fun == 'mc':
    print 'bad read_triplet_back'
    exit(1)
  else:
    # input: "figs/triplet-back-inbetween-*-%4.2f.dat" % (ff)
    filename = "figs/triplet-back-inbetween-%s-%4.2f.dat" % (fun, ff)
    data = loadtxt(filename)
  return data[:,0:4]

def read_triplet(ff, fun):
  if fun == 'mc':
    return loadtxt("figs/mc/triplet/tripletMC-%3.1f-04.20-trimmed.dat" % ff)
  elif fun == 'this-work':
    return loadtxt("figs/tripletWB-this-work-%4.2f-4.20.dat" % ff)

def read_gr(ff):
  return loadtxt("figs/gr-%04.2f.dat" % ff)

# widths given in height units (such that the figure height is 1)
# twod_width should be a constant based on figure dimensions, oned_width
# is adjustible
twod_width = 1.6
oned_width = 1.4

scale = 4
fig = figure(figsize=(scale*(twod_width + oned_width), scale))
gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[twod_width, oned_width])

xplot = subplot(gs[1])
zplot = xplot.twiny()
twod_plot = subplot(gs[0])

fig.subplots_adjust(left=0.05, right=0.975, bottom=0.15, top=0.9, wspace=0.05)

xlow = -6
xhigh = 8
xplot.set_xlim(xlow, xhigh)
zplot.set_xlim(xlow + rpath, xhigh + rpath)

xticks = [-6, -4, -2, 0]
zticks = [4, 6, 8, 10]

plotticks = xticks + [z-rpath for z in zticks]
xplot.set_xticks(plotticks)
xplot.set_xticklabels(['$%i$' %tick for tick in xticks[:-1]] + ['$0$ $0$'] +
                      ['$%i$' %tick for tick in zticks])
zplot.set_xticks([])


zplot.axvline(x=rpath, color='k')
zplot.axvline(x=3*rpath, color='k')

bracket.bracket(xplot, -0.01, -xlow/(xhigh - xlow), -.06, .06, r'$x/R$')
bracket.bracket(xplot, -xlow/(xhigh - xlow), 1.01, -.06, .06, r'$z/R$')

xmin = 1.0
xmax = 9.0
ymax = 6.0

twod_plot.set_xlim(zmin, zmax)
twod_plot.set_ylim(-rmax, rmax)

#fig.subplots_adjust(hspace=0.001)

me = 160

for name in plots:
    g3_path = read_triplet_path(ff, name)
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

    if name == 'fischer':
      # Fischer et al only predict pair distribution function in
      # contact.  We do this using "&" below which means "and".
      incontact = x**2 + (z-2*rpath)**2 < (rpath + 0.01)**2
      zplot.plot(z[incontact], g[incontact], styles.plot[name],
                 label=styles.title[name])
    else:
      myme = me
      if name == 'mc':
        myme = 10;
      xplot.plot(x[x < 0], g[x < 0], styles.plot_forward[name],
                 label=styles.title[name], markevery=myme)
      zplot.plot(z[abs(z-zcontact) > 0.0], g[abs(z-zcontact) > 0.0],
                 styles.plot_forward[name], label=styles.title[name], markevery=myme)

for name in plots:
  if name in ['this-work', 'sokolowski']:
    g3_path = read_triplet_back(ff, name)
    x = g3_path[:,3]
    z = g3_path[:,2]
    g = g3_path[:,1]
    zcontact = z.max()
    z = zcontact + (zcontact - z)
    incontact = x**2 + (z-rpath)**2 < (rpath + .01)**2

    xplot.plot(x[z==zcontact][me//2:],g[z==zcontact][me//2:], styles.plot_back[name], markevery=me)
    zplot.plot(z[z>zcontact][me//2:],g[z>zcontact][me//2:], styles.plot_back[name], markevery=me)


#xplot.set_ylabel(r'$g^{(3)}(\left< 0,0,0\right>,\left< 0,0,2.1\sigma\right>,\mathbf{r})$')
zplot.legend(loc='upper left', ncol=1).draw_frame(False)

twod_plot.set_aspect('equal')
g3mc = read_triplet(ff, 'mc')[:, center/dx:-1]
rpoints = len(g3mc[:,0])
zpoints = len(g3mc[0,:])
r = arange(0, rpoints*dx, dx)
z = arange(center, center+zpoints*dx, dx)
Z, R = meshgrid(z, r)
gmax = g3mc.max()

g3dft = read_triplet(ff, 'this-work')
zdft = loadtxt("figs/triplet-z.dat")
xdft = loadtxt("figs/triplet-x.dat")

levels = linspace(0, gmax, gmax*100)
gr = read_gr(ff)
ginf = interp(2*rpath, gr[:,0], gr[:,1]/ff)
xlo = 0.5*ginf/gmax
xhi = 1.5*ginf/gmax
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


CS = twod_plot.pcolormesh(Z, R, g3mc, vmax=gmax, vmin=0, cmap=cmap)
twod_plot.pcolormesh(-(Z-2*center), R, g3mc, vmax=gmax, vmin=0, cmap=cmap)
twod_plot.pcolormesh(zdft, -xdft, g3dft, vmax=gmax, vmin=0, cmap=cmap)
plot([zmin,zmax], [0,0], 'k-', linewidth=2)

text(-3.7, -3.9, 'this work', path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])
text(-3.7, 3.5, 'Monte Carlo', path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])

sphere0 = Circle((0, 0), 1, color='slategray')
sphere1 = Circle((2*rpath, 0), 1, color='slategray')
twod_plot.add_artist(sphere0)
twod_plot.add_artist(sphere1)

myticks = arange(0, floor(2.0*gmax)/2 + 0.1, 0.5)
colorbar(CS, extend='neither', ticks=myticks)
twod_plot.set_ylabel('$x/R$');
twod_plot.set_xlabel('$z/R$');

# Here we plot the paths on the 2d plot.  The mc plot should align
# with the dft one.
g3_path = read_triplet_path(ff, 'mc')
xmc = g3_path[:,3]
zmc = g3_path[:,2]
g3_path = read_triplet_path(ff, 'this-work')
xdft = g3_path[:,3]
zdft = g3_path[:,2]
g3_back = read_triplet_back(ff, 'this-work')
xback = g3_back[:,3]
zback = g3_back[:,2]
plot(zmc,xmc, 'w-', linewidth=3)
#plot(zdft,-xdft, 'w-', linewidth=3)
#plot(zback,-xback, 'w-', linewidth=3)
#plot(zdft,xdft, 'w-', linewidth=3)
plot(zmc,xmc, styles.color['mc']+'--', linewidth=3)
#plot(zdft,-xdft, colors[plots.index('this-work')]+'--', linewidth=3)
#plot(zback[zback<zback.max()],-xback[zback<zback.max()],
#     colors[plots.index('this-work')]+'--', linewidth=3)

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
texteff = [matplotlib.patheffects.withStroke(linewidth=2, foreground="w")]
arroweff = [matplotlib.patheffects.withStroke(linewidth=3, foreground="w")]
annotate('A', xy=(Az, Ax), xytext=(1.2,-3.5),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
annotate('B', xy=(Bz, Bx), xytext=(1.5,1.8),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
annotate('C', xy=(Cz, Cx), xytext=(rpath*2,3.0),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
annotate('D', xy=(Dz, Dx), xytext=(6.7,0.5),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
annotate('E', xy=(Ez, Ex), xytext=(8.0,0.5),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))

# Annotations on 1d plot
xplot.annotate('A', xy=(Ax, g3pathfunction_x(Ax)),
               xytext=(Ax-0.5,0.7),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
xplot.annotate('B', xy=(Bx,g3pathfunction_x(Bx)),
               xytext=(Bx+1, g3pathfunction_x(Bx)),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('C', xy=(Cz,g3pathfunction_z(Cz)),
               xytext=(Cz,g3pathfunction_z(Cz)-0.5),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('D', xy=(Dz,g3pathfunction_z(Dz)),
               xytext=(Dz+1,g3pathfunction_z(Dz)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('E', xy=(Ez,g3pathfunction_z(Ez)),
               xytext=(Ez-0.5,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))

ylim = xplot.get_ylim()
xplot.set_ylim(0, ylim[1])

plot_labels=['c)', 'd)']
# add figure labels
twod_plot.text(twod_plot.get_xlim()[0]-1, twod_plot.get_ylim()[1], plot_labels[0])
xplot.text(-7.5, ylim[1], plot_labels[1])

title = r'$g^{(3)}(\left< 0,0,0\right>,\left< 0,0,2.1\sigma\right>,\mathbf{r})$'
twod_plot.set_title(title +  ' at $\eta = %g$' % ff)
xplot.set_ylabel(title)
#fig.tight_layout(rect=[0, .03, 1, 1])
savefig("figs/triplet-correlation-pretty-inbetween-%d.pdf" % (int(ff*10)))
show()

