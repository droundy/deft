#!/usr/bin/python

from __future__ import division
from pylab import *
import matplotlib.pyplot as plt
import sys
if len(sys.argv) != 4:
    print("Usage:  " + sys.argv[0] + "-eta" + "-delta_r" + "-dv")
    exit(1)
eta = sys.argv[1]
eta = float(eta)
delta_r = sys.argv[2]
delta_r = float(delta_r)
dv = sys.argv[3]
dv = float(dv)
#pdffilename = sys.argv[1]
plt.title('$da/dz$')
plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])

def read_a1_mc():
  filename = "mc/wallsMC-pair-%02.1f-a1.dat" % eta
  print 'Using', filename
  data = loadtxt(filename)
  minr = delta_r - 3*dv
  maxr = delta_r + 3*dv
  row_min = floor(minr/0.1 + .05)
  row_max = floor(maxr/0.1 + .05)
  return sum(data[row_min:row_max], axis=0)

def read_da_dz(version):
  filename = "walls_daWB-%s-%04.2f-%04.2f-%05.3f.dat" % (version,eta,delta_r,dv) #0.%d0
  print 'Using', filename
  data = loadtxt(filename)
  z0 = data[:,0]
  da_dz = data[:,1]
  return z0, da_dz

versions = ["fischer","mc","nA","simple"]
for version in versions:
  z0, da_dz = read_da_dz(version)
  plt.plot(z0, da_dz)
plt.plot(arange(3,23,.2), read_a1_mc(), 'k')
plt.xlim([2,9])
plt.legend(versions, loc='upper right')

xlim(3, 10)

plt.show()

