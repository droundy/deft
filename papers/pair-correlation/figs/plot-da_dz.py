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

plt.figure(figsize=(6,4))

mc_z0, mc_da_dz = read_a1_mc();
plt.plot(mc_z0[::10], mc_da_dz[::10], styles.plot['mc'], label=styles.title['mc'])

for version in versions:
  z0, da_dz = read_da_dz(version)
  plt.plot(z0, da_dz, styles.plot[version], label=styles.title[version])


plt.xlim([-.5,6])
plt.legend(loc='upper right').draw_frame(False)

plt.xlabel('$z/R$')
plt.ylabel(r'$da_1/dz$')
plt.yticks([])
#plt.xticks([0,1,2,3])

if 2 < delta_r and delta_r < 2.01:
  plt.title('$\Phi(r) = \delta(\sigma - r + \delta)$     $\eta = %g$' % eta)
else:
  plt.title('$\Phi(r) = \delta(%g\sigma - r)$     $\eta = %g$' % (delta_r/2, eta))

plt.tight_layout()

savefig("figs/dadz-%d-%d.pdf" % (int(eta*10), int(delta_r)))

plt.show()
