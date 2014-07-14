#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if len(sys.argv) < 3 or sys.argv[2] != "show":
  matplotlib.use('Agg')
from pylab import *
import scipy.ndimage
import matplotlib.patheffects
import os.path
import math
import bracket # our handy bracket function
import styles # our preferred line styles

from matplotlib import rc

#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('font', family='serif', serif ='Computer Modern')
rc('text', usetex=True)

from matplotlib.colors import NoNorm

dx = 0.1
############################
rpath = 2.005

able_to_read_file = True

# Set the max parameters for plotting.
zmax = 6
zmin = -1
rmax = 4.1
############################

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])
#arg ff = [0.1, 0.2, 0.3]

plot_labels = ['a)', 'b)'] if ff == 0.1 else ['c)', 'd)']

def read_walls_path(ff, fun):
  if fun == 'mc':
    # input:  "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % (ff)
    filename = "figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" % ff
  elif fun == 'sphere-dft':
    filename = "figs/wallsWB-with-sphere-path-%1.2f.dat" % 0.3 # ff FIXME others don't exist in repo yet
  else:
    # input: "figs/walls/wallsWB-path-this-work-short-pair-%1.2f-0.005.dat" %(ff)
    # input: "figs/walls/wallsWB-path-fischer-pair-%1.2f-0.005.dat" %(ff)
    # input: "figs/walls/wallsWB-path-sokolowski-pair-%1.2f-0.005.dat" %(ff)
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-0.005.dat" %(fun, ff)
  data = loadtxt(filename)
  if fun == 'mc':
    data[:,0]-=4.995
  else:
    data[:,2]-=3
  return data[:,0:4]

def read_walls_mc(ff):
  # The 0.01 is really 0.005, but is only saved to 2 digits of precision
  return loadtxt("figs/mc/wallsMC-pair-%1.1f-0.01-trimmed.dat" % ff)

def read_walls_dft(ff, fun):
  if fun == 'sphere-dft':
    filename = "figs/wallsWB-with-sphere-%1.2f-trimmed.dat" % 0.3 # ff FIXME others don't exist in repo yet
  else:
    # input: "figs/walls/wallsWB-this-work-short-pair-%1.2f-0.05.dat" %(ff)
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-0.05.dat" %(fun, ff)
  return loadtxt(filename)

ymax = 4

# widths given in height units (such that the figure height is 1)
# twod_width should be a constant based on figure dimensions, oned_width
# is adjustible
twod_width = 1.0
oned_width = 2.0

scale = 4
fig = figure(figsize=(scale*(twod_width + oned_width), scale))
gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[twod_width, oned_width])

# xplot = fig.add_subplot(1,2,2)
# zplot = xplot.twiny()
# twod_plot = fig.add_subplot(1,2,1)

xplot = subplot(gs[1])
zplot = xplot.twiny()
twod_plot = subplot(gs[0])

fig.subplots_adjust(left=0.05, right=0.975, bottom=0.15, top=0.9, wspace=0.1)

def z_to_x(z):
  return 2 - z
def x_to_z(x):
  return 2 - x

zmax_lineplot = 6.
xmax_lineplot = 4.
xplot.set_xlim(zmax_lineplot, -xmax_lineplot)
xplot.set_xticks([6, 4, 2, 0, -2, -4])
xplot.set_xticklabels(["$6$", "$4$", "$2~~0$", "$2$", "$4$", "$6$"])
zplot.set_xlim(x_to_z(xplot.get_xlim()[0]), x_to_z(xplot.get_xlim()[1]))
zplot.set_xticks([])
#xplot.set_ylim(0)

bracket.bracket(xplot, -.01, xmax_lineplot/(xmax_lineplot+zmax_lineplot), -.06, .06, r'$x/R$')
bracket.bracket(zplot, xmax_lineplot/(xmax_lineplot+zmax_lineplot), 1.01, -.06, .06, r'$z/R$')

twod_plot.set_xlim(-1, 6)
twod_plot.set_ylim(-ymax, ymax)
sphere = Circle((0, 0), 1, color='slategray')
twod_plot.add_artist(sphere)

#fig.subplots_adjust(hspace=0.001)


twod_plot.set_aspect('equal')
g2mc = read_walls_mc(ff)
rbins = round(rmax/dx)
zposbins = round(zmax/dx)
znegbins = round(-zmin/dx)
zbins = zposbins + znegbins
g22 = zeros((rbins, zbins))
g22[:, znegbins:zbins] = g2mc[:rbins,:zposbins]
#g22[:rbins/2, znegbins:zbins] = flipud(g2mc[:rbins/2,:zposbins])
g2mc = g22
gmax = g2mc.max()
dx = 0.1

r = arange(0, rmax, dx)
z = arange(zmin, zmax, dx)
Z, R = meshgrid(z, r)

g2dft = read_walls_dft(ff, 'this-work-short')
zdft = loadtxt("figs/walls/z.dat")
xdft = loadtxt("figs/walls/x.dat")

# get rid of long range:
g2dft[xdft**2 + zdft**2 > styles.short_range**2] = 1.0
r = styles.short_range
phi = linspace(-pi/2, 0, 100)
twod_plot.plot(r*cos(phi), r*sin(phi), '-k')

levels = linspace(0, gmax, gmax*100)
xlo = 0.25/gmax
xhi = 1.25/gmax
if xhi > (1+2.0/gmax)/3:
  xhi = (1+2.0/gmax)/3
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

zextra, xextra = meshgrid(arange(zmin, -zmin, -zmin/2), arange(-rmax,rmax, rmax/2))
twod_plot.contourf(zextra, xextra, zeros_like(zextra), levels=[-1,1], colors=['k','k'])
CS = twod_plot.pcolormesh(Z, R, g2mc, vmax=gmax, vmin=0, cmap=cmap)
twod_plot.pcolormesh(zdft, -xdft, g2dft, vmax=gmax, vmin=0, cmap=cmap)

myticks = arange(0, floor(2.0*gmax)/2 + 0.1, 0.5)
fig.colorbar(CS, extend='neither', ticks=myticks)
twod_plot.set_ylabel('$x/R$');
twod_plot.set_xlabel('$z/R$');

# takes two arrays, and averages points so that a plot of x vs y
# will have points separated by a distance dpath
# returns (x, y)
def avg_points(x, y, dpath):
  new_y = array([])
  new_x = array([])
  old_i = 0
  for i in xrange(1, len(x)):
    dist = sqrt((x[i] - x[old_i])**2 + (y[i] - y[old_i])**2)
    if dist >= dpath or i == len(x) - 1:
      avg_x = average(x[old_i:i])
      avg_y = average(y[old_i:i])

      new_x = append(new_x, avg_x)
      new_y = append(new_y, avg_y)
      old_i = i
  return (new_x, new_y)


for name in styles.plots:
    g2_path = read_walls_path(ff, name)
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
    xcontact = 2.0 if name == 'mc' else 2.0051
    incontact = (x<xcontact) & (z<2)

    g_x = g[z==zcontact]
    x_x = x[z==zcontact]

    g_c = g[incontact]
    z_c = z[incontact]

    g_z = g[(x<xcontact) & (z>2.005)]
    z_z = z[(x<xcontact) & (z>2.005)]

    if 'short' in name:
        g_x = g_x[x_x <= styles.short_range]
        x_x = x_x[x_x <= styles.short_range]

        g_c = g_c[z_c <= styles.short_range]
        z_c = z_c[z_c <= styles.short_range]

        g_z = g_z[z_z <= styles.short_range]
        z_z = z_z[z_z <= styles.short_range]

    if name == 'mc':
      # do point averaging, so that points are fixed path distance apart
      dpath = 0.2
      x_x, g_x = avg_points(x_x, g_x, dpath)
      z_c, g_c = avg_points(z_c, g_c, dpath)
      z_z, g_z = avg_points(z_z, g_z, dpath)

    zplot.plot(z_c, g_c, styles.plot[name], label=styles.title[name])
    if name != 'fischer':
      # Fischer et al only predict pair distribution function in contact
      xplot.plot(x_x, g_x, styles.plot[name], label=styles.title[name])
      zplot.plot(z_z, g_z, styles.plot[name])

    # insert zoomed subplot for eta = 0.1 only
    if ff == 0.1:
      suba = axes([.65, .23, .3, .35])
      suba.plot(z_c, g_c, styles.plot[name], label=styles.title[name])
      suba.set_yticks([1, 1.1, 1.2, 1.3, 1.4])
      sub_ylim = (1, suba.get_ylim()[1])
      suba.set_ylim(sub_ylim)
      sub_xlim = suba.get_xlim()

      for i in suba.spines.itervalues():
        i.set_linewidth(2)
      if name == 'this-work-short': # only want to draw rectangle once
        zplot.add_patch(Rectangle((sub_xlim[0], sub_ylim[0]),
                                  sub_xlim[1]-sub_xlim[0],
                                  sub_ylim[1]-sub_ylim[0], facecolor='none',
                                  linewidth=2))



zplot.axvline(x=2, color='k')
zplot.axvline(x=0, color='k')


legendloc = 'lower left' if ff < 0.2 else 'upper left'
zplot.legend(loc=legendloc, ncol=1).draw_frame(False)


twod_plot.text(2.1, -3.8, styles.title['this-work-short'],
     path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])
twod_plot.text(2.1, 3.5, styles.title['mc'],
     path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])

xs = [0, 0]
ys = [ymax, rpath]
dtheta = pi/80
for theta in arange(0, pi/2 + dtheta/2, dtheta):
    xs.append(rpath*sin(theta))
    ys.append(rpath*cos(theta))
xs.append(2*ymax)

ys.append(0)
twod_plot.plot(xs, ys, 'w-', linewidth=2)
twod_plot.plot(xs, ys, 'k--', linewidth=2)


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

zplot.set_xticks([x_to_z(Ax), Bz, Cz, Dz, Ez])
zplot.set_xticklabels(["$A$", "$B$", "$C$", "$D$", "$E$"])

hw = 4 # headwidth of arrows

g2nice = read_walls_path(ff, 'this-work-short')
def g2pathfunction_x(x):
    return interp(x, flipud(g2nice[:,3]), flipud(g2nice[:,1]))
def g2pathfunction_z(z):
    return interp(z, g2nice[:,2], g2nice[:,1])

# Annotations on 2d plot
twod_plot.annotate('$A$', xy=(Az,Ax), xytext=(1,3), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
twod_plot.annotate('$B$', xy=(Bz,Bx), xytext=(1,2.5), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
twod_plot.annotate('$C$', xy=(Cz,Cx), xytext=(2.3,2.0), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
twod_plot.annotate('$D$', xy=(Dz,Dx), xytext=(3,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
twod_plot.annotate('$E$', xy=(Ez,Ex), xytext=(5,1), arrowprops=dict(shrink=0.01, width=1, headwidth=hw))

# Annotations on 1d plot
# xplot.annotate('$A$', xy=(Ax, g2pathfunction_x(Ax)),
#                xytext=(Ax-0.2, g2pathfunction_x(Ax) + ff),
#                arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
# xplot.annotate('$B$', xy=(Bx,g2pathfunction_x(Bx)),
#                xytext=(Bx+.8, g2pathfunction_x(Bx)-ff/3),
#                arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
# zplot.annotate('$C$', xy=(Cz,g2pathfunction_z(Cz)),
#                xytext=(Cz,g2pathfunction_z(Cz)-3*ff+.2),
#                arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
# zplot.annotate('$D$', xy=(Dz,g2pathfunction_z(Dz)),
#                xytext=(Dz+.7,g2pathfunction_z(Dz)+ff/2),
#                arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
# zplot.annotate('$E$', xy=(Ez,g2pathfunction_z(Ez)),
#                xytext=(Ez+0.5,g2pathfunction_z(Ez)-ff),
#                arrowprops=dict(shrink=0.01, width=1, headwidth=hw))

ylim = xplot.get_ylim()
xplot.set_ylim(0, ylim[1])


title = r'$g^{(2)}(\left< 0,0,0\right>, \left<x, 0, z\right>)$'
twod_plot.set_title(title +  r' at $\eta = %g$' % ff)
xplot.set_ylabel(title)

# add figure labels
twod_plot.text(-2.0, ymax, plot_labels[0])
xplot.text(6.7, ylim[1], plot_labels[1])

#fig.tight_layout(rect=[0, .03, 1, 1])
savefig("figs/pair-correlation-pretty-%d.pdf" % (int(ff*10)))
show()

