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


def read_mc(ff, poly, N):
  data = loadtxt("figs/mc/polyhedraMC-walls-%4.2f-density-%s-%i.dat" %(ff, poly, N))
  return data[:,0], data[:,1]

z, density = read_mc(ff, poly, N)
length = z[-1] + z[0]
print sum(density)/len(density)*length**3
plot(z, density, label="%s, %i" %(poly, N))

legend(loc='best')

title("density")
savefig("figs/density-%4.2f-%s.pdf" %(ff, poly))
show()
