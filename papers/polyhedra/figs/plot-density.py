#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os

if len(sys.argv) < 4:
  print("Use: %s ff shape N")
if len(sys.argv) < 5 or sys.argv[4] != "show":
  matplotlib.use('Agg')
from pylab import *

ff = float(sys.argv[1])
poly = sys.argv[2]
N = int(sys.argv[3])


def read_mc(ff, poly, N, dim):
  fname = "figs/mc/polyhedraMC-walls-%4.2f-%cdensity-%s-%i.dat" %(ff, dim, poly, N)
  print "using", fname
  if (not os.path.isfile(fname)):
    print("\n%s is not a file.\n\nPerhaps you have the wrong number of %ss?" %(fname, poly))
    exit(1)
  data = loadtxt(fname)
  return data[:,0], data[:,1]

dims = ['x', 'y', 'z']
#dims = ['z']
colors = ['r', 'b', 'k']
for dim, color in zip(dims, colors):
  coord, density = read_mc(ff, poly, N, dim)
  dw = coord[2] - coord[1]
  length = len(coord)*dw
  mid = length/2
  print "Integral: ", sum(density)/len(density)*length**3
  #plot(coord, density)
  plot(coord[coord<mid], density[coord<mid], label="$%c$" %dim, color=color, linestyle='-')
  plot(length-coord[coord>mid], density[coord>mid], color=color, linestyle='--')


if poly == 'cube':
  axvline(x = 0.5, linestyle=':')
  axvline(x = sqrt(2)/2, linestyle=':')
  axvline(x = sqrt(3)/2, linestyle=':')
legend(loc='best')

#xlim(0, 4)
title("density, %s, $\\eta = %04.2f$, $N = %i$." %(poly, ff, N))
savefig("figs/density-%4.2f-%s.pdf" %(ff, poly))
show()
