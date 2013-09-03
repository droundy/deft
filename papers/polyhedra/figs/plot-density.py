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
  fname = "figs/mc/polyhedraMC-walls-%4.2f-density-%s-%i.dat" %(ff, poly, N)
  if (not os.path.isfile(fname)):
    print("\n%s is not a file.\n\nPerhaps you have the wrong number of %ss?" %(fname, poly))
    exit(1)
  data = loadtxt(fname)
  return data[:,0], data[:,1]

z, density = read_mc(ff, poly, N)
dz = z[2] - z[1]
length = len(z)*dz
zmid = length/2
print sum(density)/len(density)*length**3
plot(z[z<zmid], density[z<zmid], label="%s, $\\eta = %04.2f$,  $N = %i$" %(poly, ff, N))
plot(length-z[z>zmid], density[z>zmid], label="%s, $\\eta = %04.2f$,  $N = %i$" %(poly, ff, N))
if poly == 'cube':
  axvline(x = 0.5, linestyle=':')
  axvline(x = sqrt(2)/2, linestyle=':')
  axvline(x = sqrt(3)/2, linestyle=':')
legend(loc='best')

#xlim(0, 4)
title("density")
savefig("figs/density-%4.2f-%s.pdf" %(ff, poly))
show()
