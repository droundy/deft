#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os
import common

if len(sys.argv) < 4:
  print("Use: %s ff shape N")
if len(sys.argv) < 5 or sys.argv[4] != "show":
  matplotlib.use('Agg')
from pylab import *

ff = float(sys.argv[1])
poly = sys.argv[2]
N = int(sys.argv[3])


e, densities = common.read_mc_density(ff, poly, N, 'walls')
length = common.read_mc_dimensions(ff, poly, N, 'walls')
print length
dims = ['x', 'y', 'z']
colors = ['r', 'b', 'k']
for i in xrange(3):
  coord = e[i]
  density = densities[i]
  mid = length[i]/2
  print("Integral %c: %g" %(dims[i], sum(density)/len(density)*length[0]*length[1]*length[2]))
  plot(coord[coord<mid], density[coord<mid], label="$%c$" %dims[i], color=colors[i], linestyle='-')
  plot(length[i]-coord[coord>mid], density[coord>mid], color=colors[i], linestyle='--')


if poly == 'cube':
  axvline(x = 0.5, linestyle=':')
  axvline(x = sqrt(2)/2, linestyle=':')
  axvline(x = sqrt(3)/2, linestyle=':')
legend(loc='best')

#xlim(0, 4)
title("density, %s, $\\eta = %04.2f$, $N = %i$." %(poly, ff, N))
savefig("figs/density-%4.2f-%s.pdf" %(ff, poly))
show()
