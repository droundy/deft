#!/usr/bin/python2
from __future__ import division
import matplotlib, sys, os, argparse
import read, arguments

parser = argparse.ArgumentParser(
  description='Plot density of polyhedra.', parents = [arguments.parser])

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
  N = read.get_N("figs/mc/%s-%4.2f-g-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N

if args.hide:
  matplotlib.use('Agg')
from pylab import *


names, data = read.read_mc_g(ff, polyhedron, N, celltype)
colors = ['k', 'b', 'c', 'r', 'g', 'g']
xmax = max(data[:,0])
xlim(0,xmax)

for i in xrange(1, len(names)):
  plot(data[:,0], data[:,i], label=names[i], color=colors[i-1])

legend(loc='best')
title("distribution functions, %s, %s, $\\eta = %04.2f$, $N = %i$." %(celltype, polyhedron, ff, N))
savefig("figs/%s-g-%4.2f-%s.pdf" %(celltype, ff, polyhedron))
show()
