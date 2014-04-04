#!/usr/bin/python2
from __future__ import division
import matplotlib, sys, os, argparse
import read

parser = argparse.ArgumentParser(description='Plot order parameter of polyhedra.')
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
args = parser.parse_args()


ff = args.ff
polyhedron = args.shape

if args.periodic:
  celltype = 'periodic'
else:
  celltype = 'walls'

if args.N == 0:
  N = read.get_N("figs/mc/%s-%4.2f-order-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N

if args.hide:
  matplotlib.use('Agg')
from pylab import *

order_parameters = read.read_mc_order(ff, polyhedron, N, celltype)

dim = read.read_mc_dimensions(ff, polyhedron, N, celltype)
dz = dim[2]/len(order_parameters[0,:])
dcostheta = 1/len(order_parameters[:,0])
print dz
print dcostheta

z = arange(0, dim[2], dz)
costheta = arange(0, 1, dcostheta)
z, costheta = meshgrid(z, costheta)

ylim(cos(arctan(sqrt(2))), 1)
print order_parameters.shape
print z.shape

highest = nanmax(order_parameters.flat)
print highest
pcolormesh(z, costheta, order_parameters, vmax=highest/10, vmin=0)
print len(z)

show()
