#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os, argparse
import common

parser = argparse.ArgumentParser(description='Plot pressure of polyhedra.')
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
  N = common.get_N("figs/mc/%s-%4.2f-pressure-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N

if args.hide:
  matplotlib.use('Agg')
from pylab import *

dV, totalmoves, pressure, dZ = common.read_mc_pressure(ff, polyhedron, N, celltype)

avgpressure = ones_like(pressure)*(sum(dZ)/dV/totalmoves[-1])

plot(totalmoves/N, avgpressure, "-", label="avg")
plot(totalmoves/N, pressure, "o")

title("$%s,$ $pressure,$ $%s,$ $\\eta = %04.2f$, $N = %i$." %(celltype, polyhedron, ff, N))
savefig("figs/%s-pressure-%4.2f-%s.pdf" %(celltype, ff, polyhedron))
show()
