#!/usr/bin/python

from __future__ import division
import matplotlib
import sys
import styles # for our style choices for these plots
if not(len(sys.argv) >= 4 and sys.argv[3] == "show"):
  matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
if len(sys.argv) < 3:
    print("Usage:  " + sys.argv[0] + " eta delta_r")
    exit(1)
eta = float(sys.argv[1])
#arg eta = [.2,.3]
delta_r = float(sys.argv[2])
#arg delta_r = [2.005, 3.005]
able_to_read_file = True

plt.title('$da/dz,$ $\eta = %g,$  $\Delta r = %g$' %(eta, delta_r))

def read_a1_mc():
  data = loadtxt("figs/mc/a1/wallsMC-a1-pair-%02.1f-%1.3f.dat" %(eta,delta_r))
  return data[:,0], data[:,1]

def read_da_dz(version):
  # input: "figs/walls.dat" % ()
  # input: "figs/walls/walls_daWB-*-%04.2f-%05.3f.dat" % (eta,delta_r)
  filename = "figs/walls/walls_daWB-%s-%04.2f-%05.3f.dat" % (version,eta,delta_r) #0.%d0
  #print 'Using', filename
  try:
    data = loadtxt(filename)
  except IOError:
    global able_to_read_file
    able_to_read_file = False
    print "File not found: ", filename
    return 0,0
  z0 = data[:,0]
  da_dz = data[:,1]
  return z0, da_dz

versions = ["this-work", 'sokolowski', "fischer"]

plt.figure(figsize=(5,3))

for version in versions:
  z0, da_dz = read_da_dz(version)
  plt.plot(z0/2, da_dz, styles.plot[version], label=styles.title[version])

mc_z0, mc_da_dz = read_a1_mc();


plt.plot(mc_z0/2, mc_da_dz, styles.plot['mc'], label=styles.title['mc'])
plt.xlim([-.5/2.,6/2.])
plt.legend(loc='best').draw_frame(False)

plt.xlabel('$z/\sigma$')
plt.ylabel(r'$\frac{da_1}{dz}$')
plt.yticks([])
plt.xticks([0,1,2,3])
plt.tight_layout()

#plotname = "figs/dadz-" + str(int(eta*10)) + "-" + str(int(delta_r)) + ".pdf"
savefig("figs/dadz-%d-%d.pdf" % (int(eta*10), int(delta_r)))

plt.show()
