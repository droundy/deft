#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if len(sys.argv) < 3 or sys.argv[2] != "show":
  matplotlib.use('Agg')
from pylab import *
import scipy.ndimage
import os.path
import math
import matplotlib.patheffects
import bracket # our handy bracket function
import styles # our preferred line styles
from matplotlib import rc
rc('text', usetex=True)

from matplotlib.colors import NoNorm

# these are the things to set
plots = ['mc', 'this-work', 'this-work-mc', 'sokolowski', 'fischer']#, 'sphere-dft'] # , 'gloor'
dx = 0.1
sigma = 2.0
contact_delta = 0.1
############################
rpath = sigma + contact_delta
center = rpath/2


able_to_read_file = True

# Set the max parameters for plotting.
zmax = 7
zmin = -3
rmax = 4
############################

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])
#arg ff = [0.1, 0.2, 0.3, 0.4]

def read_triplet_path(ff, fun):
  if fun == 'mc':
    # input:  "figs/mc/triplet/tripletMC-%03.1f-path-trimmed.dat" % (ff)
    filename = "figs/mc/triplet/tripletMC-%03.1f-path-trimmed.dat" % (ff)
  else:
    # input: "figs/tripletWB-path-*-%1.2f.dat" %(ff)
    filename = "figs/tripletWB-path-%s-%1.2f.dat" %(fun, ff)
  if (os.path.isfile(filename) == False):
    # MC data is in repo, but dft isn't, so just use that for now so it will build.
    # input:  "figs/mc/triplet/tripletMC-%03.1f-path-trimmed.dat" % (ff)
    filename = "figs/mc/triplet/tripletMC-%03.1f-path-trimmed.dat" % (ff)

  data = loadtxt(filename)
  if fun == 'mc':
    data[:,0]-=4.995
  #data[:,0:4][:,1] = data[:,0:4][:,1]*0 + 1 # fixme remove
  return data[:,0:4]

def read_triplet_back(ff, fun):
  if fun == 'mc':
    print 'bad read_triplet_back'
    exit(1)
  else:
    # input: "figs/triplet-back-contact-*-%4.2f.dat" % (ff)
    filename = "figs/triplet-back-contact-%s-%4.2f.dat" % (fun, ff)
    data = loadtxt(filename)
  #data[:,0:4][:,1] = data[:,0:4][:,1]*0 + 1 # fixme remove
  return data[:,0:4]

def read_triplet(ff, fun):
  if fun == 'mc':
    return loadtxt("figs/mc/triplet/tripletMC-%3.1f-02.10-trimmed.dat" % ff)
  elif fun == 'this-work':
    return loadtxt("figs/tripletWB-this-work-%4.2f-2.10.dat" % ff)

def read_gr(ff):
  return loadtxt("figs/gr-%04.2f.dat" % ff)

# widths given in height units (such that the figure height is 1)
# twod_width should be a constant based on figure dimensions, oned_width
# is adjustible
twod_width = 1.3
oned_width = 1.7

scale = 4
fig = figure(figsize=(scale*(twod_width + oned_width), scale))
gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[twod_width, oned_width])

xplot = subplot(gs[1])
zplot = xplot.twiny()
twod_plot = subplot(gs[0])

fig.subplots_adjust(left=0.05, right=0.975, bottom=0.15, top=0.9, wspace=0.05)

xmin = rpath/2*sqrt(3)
xlow = 6
xhigh = -6+xmin
xplot.set_xlim(xlow, xhigh)
zplot.set_xlim(rpath/2 - xlow + xmin, rpath/2 - xhigh + xmin)

xticks = [6, 4, 2]
zticks = [1, 3, 5, 7]

plotticks = xticks + [rpath/2 - z + xmin for z in zticks]
xplot.set_xticks(xticks)
zplot.set_xticks(zticks)
xplot.set_xticks(plotticks)
xplot.set_xticklabels(['$%i$' %tick for tick in xticks] +
                       ['$%i$' %tick for tick in zticks])
zplot.set_xticks([])

bracket.bracket(xplot, -.01, (xlow - 2)/(xlow - xhigh) + 0.007, -.06, .06, r'$x/R$')
bracket.bracket(xplot, (xlow - xmin)/(xlow - xhigh) - 0.013, 1.01, -.06, .06, r'$z/R$')

twod_plot.set_xlim(zmin, zmax)
twod_plot.set_ylim(-rmax, rmax)

#fig.subplots_adjust(hspace=0.001)

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
ginf = interp(rpath, gr[:,0], gr[:,1]/ff)
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
twod_plot.plot([zmin,zmax], [0,0], 'k-', linewidth=2)

twod_plot.text(-2.7, -3.9, 'this work', path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])
twod_plot.text(-2.7, 3.5, 'Monte Carlo', path_effects=[matplotlib.patheffects.withStroke(linewidth=2, foreground="w")])

sphere0 = Circle((0, 0), 1, color='slategray')
sphere1 = Circle((rpath, 0), 1, color='slategray')
twod_plot.add_artist(sphere0)
twod_plot.add_artist(sphere1)

myticks = arange(0, floor(2.0*gmax)/2 + 0.1, 0.5)
fig.colorbar(CS, extend='neither') # , ticks=myticks)
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
  incontact = x**2 + (z-rpath)**2 < (rpath + .01)**2

  g_x = g[z==zcontact]
  x_x = x[z==zcontact]

  g_c = g[incontact][:-1]
  z_c = z[incontact][:-1]

  g_z = g[z>rpath*2]
  z_z = z[z>rpath*2]

  if name == 'mc':
    # do point averaging, so that points are fixed path distance apart
    dpath = 0.3
    x_x, g_x = avg_points(x_x, g_x, dpath)
    z_c, g_c = avg_points(z_c, g_c, dpath)
    z_z, g_z = avg_points(z_z, g_z, dpath)

  zplot.plot(z_c, g_c, styles.plot[name], label=styles.title[name])
  if name != 'fisher':
    # Fischer et al only predict pair distribution function in contact
    xplot.plot(x_x, g_x, styles.plot[name], label=styles.title[name])
    zplot.plot(z_z, g_z, styles.plot[name], label=styles.title[name])

  # insert zoomed-in subplot
  suba = axes([.73, .59, .23, .28])
  suba.plot(z_c, g_c, styles.plot[name], label=styles.title[name])
  sub_ylim = (3.5, 5)
  suba.set_yticks(arange(sub_ylim[0], sub_ylim[1]+0.25, 0.5))
  suba.set_ylim(sub_ylim)
  sub_xlim = (1.5, z_c.max())
  suba.set_xlim(sub_xlim)

  zplot.add_patch(Rectangle((sub_xlim[0], sub_ylim[0]),
                            sub_xlim[1]-sub_xlim[0], sub_ylim[1]-sub_ylim[0], facecolor='none',
                            linewidth=2))
  for i in suba.spines.itervalues():
    i.set_linewidth(2)

xplot.axvline(x=xmin, color='k')
zplot.axvline(x=2*rpath, color='k')


def get_closest(x, y, xref, space):
  xnew = zeros_like(xref)
  ynew = zeros_like(xref)
  for i in xrange(len(xnew)):
    cond = abs((x - xref[i]) - space) == abs((x - xref[i]) - space).min()
    xnew[i] = x[cond][0]
    ynew[i] = y[cond][0]
  return (xnew, ynew)

# add triangles to show direction of plots
spacing = 2

# Make reference points
x_ref = arange(xmin + 0.1, xlow, spacing)
z_ref = arange(rpath/2 + 0.1 + spacing/2, 6+spacing, spacing)

for name in ['this-work', 'this-work-mc', 'sokolowski']:
  # forward arrows
  g3_path = read_triplet_path(ff, name)
  x = g3_path[:,3]
  z = g3_path[:,2]
  g = g3_path[:,1]
  zcontact = z.min()

  g_x = g[z==zcontact]
  x_x = x[z==zcontact]
  g_z = g[z>zcontact][::-1]
  z_z = z[z>zcontact][::-1]

  start = styles.start[name]*spacing
  x_x, g_x = get_closest(x_x, g_x, x_ref, start)
  xplot.plot(x_x, g_x, styles.plot_forward[name], mec='none')
  #if name == 'this-work': print z_z
  z_z, g_z = get_closest(z_z, g_z, z_ref, start)
  #print name, z_z
  zplot.plot(z_z, g_z, styles.plot_forward[name], mec='none')
  suba.plot(z_z, g_z, styles.plot_forward[name], mec='none')

for name in ['this-work', 'sokolowski']:
  # backward arrows
  g3_path = read_triplet_back(ff, name)
  x = g3_path[:,3]
  z = g3_path[:,2]
  g = g3_path[:,1]
  zcontact = z.max()
  z = zcontact + (zcontact - z)

  g_x = g[z==zcontact]
  x_x = x[z==zcontact]
  g_z = g[z>zcontact]
  z_z = z[z>zcontact]

  start = styles.start[name]*spacing + 0.5*spacing

  xplot.plot(x_x, g_x, styles.plot[name])
  x_x, g_x = get_closest(x_x, g_x, x_ref, start)
  xplot.plot(x_x, g_x, styles.plot_back[name], mec='none')

  start = styles.start[name]*spacing - 0.5*spacing

  zplot.plot(z_z, g_z, styles.plot[name])
  suba.plot(z_z, g_z, styles.plot[name])
  z_z, g_z = get_closest(z_z, g_z, z_ref, start)
  zplot.plot(z_z, g_z, styles.plot_back[name], mec='none')
  suba.plot(z_z, g_z, styles.plot_back[name], mec='none')






#xplot.set_ylabel(r'$g^{(3)}(\left< 0,0,0\right>,\left< 0,0,1.1\sigma\right>,\mathbf{r})$')
xplot.legend(loc='upper left', ncol=1).draw_frame(False)


Ax = 3.9
Az = rpath/2
Bx = xmin
Bz = rpath/2
Cx = rpath*(sqrt(3)/2)
Cz = rpath*1.5
Dx = 0
Dz = rpath*2
Ex = 0
Ez = 6

hw = 4 # headwidth of arrows

g3nice = read_triplet_path(ff, 'this-work')
znice = g3nice[:,2]
xnice = g3nice[:,3]
gnice = g3nice[:,1]
znicecontact = znice.min()
def g3pathfunction_x(x):
    return interp(x, xnice[znice == znicecontact], gnice[znice == znicecontact])
def g3pathfunction_z(z):
    return interp(z, flipud(znice), flipud(gnice))

# Annotations on 2d plot
texteff = [matplotlib.patheffects.withStroke(linewidth=2, foreground="w")]
arroweff = [matplotlib.patheffects.withStroke(linewidth=3, foreground="w")]
twod_plot.annotate('A', xy=(Az, Ax), xytext=(0.2,3),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
twod_plot.annotate('B', xy=(Bz, Bx), xytext=(rpath,3),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
twod_plot.annotate('C', xy=(Cz, Cx), xytext=(4.3,2.5),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
twod_plot.annotate('D', xy=(Dz, Dx), xytext=(5,1),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))
twod_plot.annotate('E', xy=(Ez, Ex), xytext=(6,1),
         path_effects=texteff,
         arrowprops=dict(shrink=0.01, width=1, headwidth=hw, path_effects=arroweff))

# Annotations on 1d plot
xplot.annotate('A', xy=(Ax, g3pathfunction_x(Ax)),
               xytext=(Ax+1,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
xplot.annotate('B', xy=(Bx,g3pathfunction_x(Bx)),
               xytext=(Bx+1, g3pathfunction_x(Bx)-1.0),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('C', xy=(Cz,g3pathfunction_z(Cz)),
               xytext=(Cz,g3pathfunction_z(Cz)-2.0),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('D', xy=(Dz,g3pathfunction_z(Dz)),
               xytext=(Dz+1,g3pathfunction_z(Dz)-0.2),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))
zplot.annotate('E', xy=(Ez,g3pathfunction_z(Ez)),
               xytext=(Ez+0.5,1.3),
               arrowprops=dict(shrink=0.01, width=1, headwidth=hw))

ylim = xplot.get_ylim()
xplot.set_ylim(0, ylim[1])

plot_labels=['a)', 'b)']
# add figure labels
twod_plot.text(twod_plot.get_xlim()[0]-1, twod_plot.get_ylim()[1], plot_labels[0])
xplot.text(6.8, ylim[1], plot_labels[1])

title = r'$g^{(3)}(\left< 0,0,0\right>,\left< 0,0,\sigma\right>,\mathbf{r})$'
twod_plot.set_title(title + ' at $\eta = %g$' % ff)
xplot.set_ylabel(title)
#fig.tight_layout(rect=[0, .03, 1, 1])
savefig("figs/triplet-correlation-pretty-contact-%d.pdf" % (int(ff*10)))
show()

