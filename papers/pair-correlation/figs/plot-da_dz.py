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

plt.xlim([2,9])
plt.legend(versions, loc='upper left')

plt.show()

