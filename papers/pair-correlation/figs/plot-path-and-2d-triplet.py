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

# these are the things to set
colors = ['k', 'b', 'g', 'r']
plots = ['mc']#, 'this-work', 'fischer'] # , 'gloor'
titles = ['Monte Carlo']#, 'this work', 'Fischer et al'] # , 'gloor'
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

def read_walls_path(ff,z0,fun):
  if fun == 'mc':
    filename = "figs/mc/triplet/wallsMC-triplet-%02.1f-path.dat" % ff
  else:
    filename = "figs/walls/wallsWB-path-%s-triplet-%1.2f-%1.2f.dat" %(fun, ff, z0)
  if (os.path.isfile(filename) == False):
    # Just use walls data if we do not have the MC (need to be careful!)
    filename = "figs/walls/wallsWB-path-%s-triplet-%1.2f-%1.2f.dat" %('this-work', ff, z0)
  data=numpy.loadtxt(filename)
  #if fun == 'mc':
  #  data[:,0]-=4.995
  return data[:,0:2]

def read_walls(ff, z0, fun):
  if fun == 'mc':
    filename = "figs/mc/triplet/wallsMC-triplet-%1.1f-%1.2f.dat" % (ff, z0)
  else:
    filename = "figs/walls/wallsWB-%s-triplet-%1.2f-%1.2f.dat" %(fun, ff, z0)
  #print 'Using', filename
  if (os.path.isfile(filename) == False):
    # Just use walls data if we do not have the MC (need to be careful!)
    filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %('this-work', ff, z0)
    pylab.title("Using this work instead of MC!")
  data = numpy.loadtxt(filename)
  return data


data = read_walls_path(ff, z0, plots[0])
print data
pylab.plot(data[:,0], data[:,1])

pylab.show()
