#!/usr/bin/python

from __future__ import division
import matplotlib, sys
if len(sys.argv) < 3 or sys.argv[2] != "show":
  matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab, numpy, scipy.ndimage
import os.path
import math

from matplotlib.colors import NoNorm

at_wall = False

# these are the things to set
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
#arg ff = [0.1, 0.2, 0.3, 0.4, 0.5]

def read_walls(ff, z0):
  # input: "figs/mc/wallsMC-pair-%1.1f-*.dat" % (ff)
  filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
  data = numpy.loadtxt(filename)
  return data

def read_g2_and_prepare_plot(zmin, zmax, rmax, centered = True):
  g2mc = read_walls(ff, z0)
  rbins = 2*round(rmax/dx)
  if not centered:
    rbins = round(rmax/dx)
  zbins = round((zmax-zmin)/dx)
  zposbins = round(zmax/dx)
  znegbins = zbins - zposbins
  g22 = numpy.zeros((rbins, zbins) )
  if centered:
    g22[rbins/2:rbins, znegbins:zbins] = g2mc[:rbins/2,:zposbins]
    g22[:rbins/2, znegbins:zbins] = numpy.flipud(g2mc[:rbins/2,:zposbins])
    g2mc = g22
  gmax = g2mc.max()

  r = numpy.linspace(0-rmax, rmax, rbins)
  z = numpy.linspace(zmin, zmax, zbins)
  Z, R = numpy.meshgrid(z, r)

  levels = numpy.linspace(0, gmax, gmax*100)
  xlo = 0.85/gmax
  xhi = 1.15/gmax
  xhier = (1 + xhi)/2.0
  whiterange = 0.03

  cdict = {'red':   [(0.0,  0.0, 0.0),
                     (xlo,  1.0, 1.0),
                     ((1-whiterange)/gmax,  1.0, 1.0),
                     ((1+whiterange)/gmax,  1.0, 1.0),
                     (xhi,  0, 0),
                     (xhier,0.0, 0.0),
                     (1.0,  1.0, 1.0)],
           'green': [(0.0, 0.0, 0.0),
                     (xlo,  0.1, 0.1),
                     ((1-whiterange)/gmax,  1.0, 1.0),
                     ((1+whiterange)/gmax,  1.0, 1.0),
                     (xhi, 0, 0),
                     (xhier,1.0, 1.0),
                     (1.0, 1.0, 1.0)],
           'blue':  [(0.0,  0.0, 0.0),
                     (xlo,  0.1, 0.1),
                     ((1-whiterange)/gmax,  1.0, 1.0),
                     ((1+whiterange)/gmax,  1.0, 1.0),
                     (xhi,  1.0, 1.0),
                     (xhier,0.0, 0.0),
                     (1.0,  0.0, 0.0)]}
  cmap = matplotlib.colors.LinearSegmentedColormap('mine', cdict)
  return Z, R, g2mc, cmap, gmax

pylab.figure(figsize=(21,15 + 10*dx))
pylab.subplot(1,1,1).set_aspect('equal')

Z, R, g2mc, cmap, gmax = read_g2_and_prepare_plot(zmin, 21, 15.0)
rbins = Z.shape[0]
nrmin = round(rbins/2) + 10
levels = pylab.concatenate([pylab.linspace(0, .85, 10),
                            pylab.linspace(.855, 0.99, 10),
                            pylab.linspace(1.01, 1.15, 10),
                            pylab.linspace(1.2, gmax, 10)])
g2mc = (1+pylab.tanh((R+10)/2))/2*(g2mc-1) + 1
pylab.contourf(Z[:nrmin,:], R[:nrmin,:], g2mc[:nrmin,:], levels, vmax=gmax*.5, vmin=0, cmap=cmap)
pylab.axes().get_xaxis().set_visible(False)
pylab.axes().get_yaxis().set_visible(False)
pylab.title('')
pylab.subplots_adjust(left=-0.01, right=1.02, top=1.02, bottom=-0.02)
pylab.savefig("figs/pretty-%d.svg"  % (int(ff*10)))

pylab.show()

