#!/usr/bin/python
from __future__ import division
import matplotlib, sys, os

if len(sys.argv) < 3:
  print("Use: %s shape N")
if len(sys.argv) < 4 or sys.argv[3] != "show":
  matplotlib.use('Agg')

from pylab import *

poly = sys.argv[1]
N = sys.argv[2]


def read_mc(poly, N):
  # input: "figs/mc/wallsMC-density-%s-%s.dat" %(poly, N)
  data = loadtxt("figs/mc/wallsMC-density-%s-%s.dat" %(poly, N))
  return data[:,0], data[:,1]

z, density = read_mc(poly, N)

plot(z, density, label=poly+', '+N)

legend(loc='best')

title("density")
savefig("figs/density-%s-%s.pdf" %(poly, N))
show()
