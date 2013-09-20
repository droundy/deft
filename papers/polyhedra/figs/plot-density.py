#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os, argparse
import read

parser = argparse.ArgumentParser(description='Plot density of polyhedra.')
parser.add_argument('ff', metavar='ff', type=float, help='filling fraction')
parser.add_argument('-N', metavar='N', type=int,default=0,
                    help="""number of polyhedra, if not supplied then the first
                            file with the proper filling fraction will be used""")
parser.add_argument('-s', '--shape', metavar='S', default='truncated_tetrahedron',
                    choices=['cube', 'tetrahedron', 'truncated_tetrahedron'],
                    help='type of polyhedron, defaults to truncated_tetrahedron')
parser.add_argument('-p', '--periodic', action='store_true',
                    help='will use periodic cell - defaults to walls otherwise')
parser.add_argument('--hide', action='store_true',
                    help='will just save the plot and won\'t display it')
parser.add_argument('-z', action='store_true', help='only plot the z-dimension')

args = parser.parse_args()


ff = args.ff
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
if args.z:
  coord = e[2]
  density=densities[2]
  plot(coord, density, label="$%c$" %dims[2], color=colors[2], linestyle='-')
else:
  for i in xrange(3):
    coord = e[i]
    density = densities[i]
    mid = length[i]/2
    print("Integral %c: %g" %(dims[i], sum(density)/len(density)*length[0]*length[1]*length[2]))
    plot(coord[coord<mid], density[coord<mid], label="$%c$" %dims[i], color=colors[i], linestyle='-')
    plot(length[i]-coord[coord>mid], density[coord>mid], color=colors[i], linestyle='--')


if polyhedron == 'cube' and celltype == 'walls':
  axvline(x = 0.5, linestyle=':')
  axvline(x = sqrt(2)/2, linestyle=':')
  axvline(x = sqrt(3)/2, linestyle=':')
  if args.z:
    axvline(x = length[2]-0.5, linestyle=':')
    axvline(x = length[2]-sqrt(2)/2, linestyle=':')
    axvline(x = length[2]-sqrt(3)/2, linestyle=':')
if not args.z:
  legend(loc='best')

#xlim(0, 4)
title("$%s,$ $density,$ $%s,$ $\\eta = %04.2f$, $N = %i$." %(celltype, polyhedron, ff, N))
savefig("figs/%s-density-%4.2f-%s.pdf" %(celltype, ff, polyhedron))
show()
