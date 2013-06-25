#!/usr/bin/python

from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import sys
if len(sys.argv) != 3:
    print("Usage:  " + sys.argv[0] + " eta delta_r")
    exit(1)
eta = sys.argv[1]
eta = float(eta)
delta_r = sys.argv[2]
delta_r = float(delta_r)
able_to_read_file = True

plt.title('$da/dz$')
plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow'])

def read_a1_mc():
  filename = "figs/mc/wallsMC-pair-%02.1f-a1.dat" % eta
  print 'Using', filename
  try:
    data = loadtxt(filename)
  except IOError:
    global able_to_read_file
    able_to_read_file = False
    return 0
  row_min = floor(delta_r/0.1-.05)
  row_max = row_min+2
  return sum(data[row_min:row_max], axis=0)/(row_max-row_min)

def read_da_dz(version):
  filename = "figs/walls/walls_daWB-%s-%04.2f-%04.2f.dat" % (version,eta,delta_r) #0.%d0
  'Using', filename
  try:
    data = loadtxt(filename)
  except IOError:
    global able_to_read_file
    able_to_read_file = False
    return 0,0
  z0 = data[:,0]
  da_dz = data[:,1]
  return z0, da_dz

versions = ["fischer","gross","this-work"]

for version in versions:
  z0, da_dz = read_da_dz(version)
  plt.plot(z0, da_dz)

mc = read_a1_mc();


print able_to_read_file


if (able_to_read_file == False):
    plt.plot(arange(0,10,1), [0]*10, 'k')
    plt.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=20)
else:
    plt.plot(arange(3.05,3.05+len(mc)*0.1,.1), mc, 'k')
    plt.xlim([2,9])
    plt.legend(versions+['mc'], loc='upper right')

plotname = "figs/dadz-" + str(int(eta*10)) + "-" + str(int(delta_r)) + ".pdf"
savefig(plotname)

plt.show()

