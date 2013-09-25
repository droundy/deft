#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os, argparse
import read, arguments

parser = argparse.ArgumentParser(
  description='Plot density of polyhedra.', parents = [arguments.parser])
parser.add_argument(
  '--reflect', action='store_true',
  help='will reflect plots halfway across cell')

args = parser.parse_args()


ff = args.ff
if args.ratio != 0:
  polyhedron = args.shape + "_%05.2f" %args.ratio
else:
  polyhedron = args.shape

if args.periodic:
  celltype = 'periodic'
else:
  celltype = 'walls'

if args.N == 0:
  N = read.get_N("figs/mc/%s-%4.2f-density-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N

if args.hide:
  matplotlib.use('Agg')
from pylab import *

e, densities = read.read_mc_density(ff, polyhedron, N, celltype)
if e == 0:
  exit(1)
length = read.read_mc_dimensions(ff, polyhedron, N, celltype)
print 'cell shape: ', length
dims = ['x', 'y', 'z']
colors = ['r', 'b', 'k']
for i in xrange(3):
  coord = e[i]
  density = densities[i]
  if args.reflect:
    mid = length[i]/2
  else:
    mid = len(coord)
  print("Integral %c: %g" %(dims[i], sum(density)/len(density)*length[0]*length[1]*length[2]))
  plot(coord[coord<mid], density[coord<mid], label="$%c$" %dims[i], color=colors[i], linestyle='-')
  if args.reflect:
    plot(length[i]-coord[coord>mid], density[coord>mid], color=colors[i], linestyle='--')


if polyhedron == 'cube' and celltype == 'walls':
  axvline(x = 0.5, linestyle=':')
  axvline(x = sqrt(2)/2, linestyle=':')
  axvline(x = sqrt(3)/2, linestyle=':')
legend(loc='best')

#xlim(0, 4)
title("density, %s, %s, $\\eta = %04.2f$, $N = %i$." %(celltype, polyhedron, ff, N))
savefig("figs/%s-density-%4.2f-%s.pdf" %(celltype, ff, polyhedron))
show()
