#!/usr/bin/python

from __future__ import division
import matplotlib
import sys
if not 'show' in sys.argv:
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
dz = 0.01

plt.figure(figsize=(4,4))
plt.title(r'$\Phi(r) = \delta(r - %g\, \sigma)$' %(delta_r/2))

def read_a1_mc():
  data = loadtxt("../../papers/pair-correlation/figs/mc/a1/wallsMC-a1-pair-%02.1f-%1.3f.dat" %(eta,delta_r))
  #print 'Using', filename
  return data[:,0], data[:,1]

def read_da_dz(version):
  # input: "../../papers/pair-correlation/figs/walls/walls_daWB-this-work-short-%04.2f-%05.3f.dat" % (eta,delta_r)
  # input: "../../papers/pair-correlation/figs/walls/walls_daWB-sokolowski-%04.2f-%05.3f.dat" % (eta,delta_r)
  # input: "../../papers/pair-correlation/figs/walls/walls_daWB-fischer-%04.2f-%05.3f.dat" % (eta,delta_r)
  filename = "../../papers/pair-correlation/figs/walls/walls_daWB-%s-%04.2f-%05.3f.dat" % (version,eta,delta_r) #0.%d0
  data = loadtxt(filename)
  z0 = data[:,0]
  da_dz = data[:,1]
  return z0, da_dz

versions = ["this-work-short","sokolowski", "fischer"]
names = ['this work', 'Sokolowski', 'Fischer']
colors = ['b', 'r', 'g']
if delta_r > 2.1:
  fischeri = versions.index('fischer')
  versions.remove('fischer')
  names.pop(fischeri)
  colors.pop(fischeri)

mc_z0, mc_da_dz = read_a1_mc();
plt.plot(mc_z0/2, np.pi*delta_r**2*mc_da_dz, 'k.', label='MC')

for i in xrange(len(versions)):
  z0, da_dz = read_da_dz(versions[i])
  plt.plot(z0/2, np.pi*delta_r**2*da_dz, colors[i]+'-', label=names[i])

plt.xlim([-.2,1.5])
#plt.ylabel(r'$da_1/dz$')
plt.xlabel(r'$z/\sigma$')
plt.legend(loc='best').draw_frame(False)

#plotname = "figs/dadz-" + str(int(eta*10)) + "-" + str(int(delta_r)) + ".pdf"
savefig("figs/dadz-%d-%d.pdf" % (int(eta*10), int(delta_r)), transparent=True)

plt.show()
