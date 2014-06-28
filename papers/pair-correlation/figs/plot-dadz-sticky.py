#!/usr/bin/python2
from __future__ import division
import matplotlib
import sys
import styles # for our style choices for these plots
if "show" not in sys.argv:
  matplotlib.use('Agg')
from pylab import *

if len(sys.argv) < 3:
    print("Usage:  " + sys.argv[0] + " eta delta_r")
    exit(1)

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

eta = float(sys.argv[1])
#arg eta = [.2,.3]
delta_r = float(sys.argv[2])
#arg delta_r = [2.005, 3.005]

z = {}
dadz = {}

versions = ['mc', 'this-work-short', 'sokolowski', 'fischer']

data = loadtxt("figs/mc/a1/wallsMC-a1-pair-%02.1f-%1.3f.dat" %(eta,delta_r))
z['mc'], dadz['mc'] = data[:,0][::10], data[:,1][::10]

# data = loadtxt("figs/walls/walls_daWB-this-work-FIXME-%04.2f-%05.3f.dat" % (eta, delta_r))
# z['this-work'], dadz['this-work'] = data[:,0], data[:,1]

# data = loadtxt("figs/walls/walls_daWB-this-work-mc-%04.2f-%05.3f.dat" % (eta, delta_r))
# z['this-work-mc'], dadz['this-work-mc'] = data[:,0], data[:,1]

data = loadtxt("figs/walls/walls_daWB-this-work-short-%04.2f-%05.3f.dat" % (eta, delta_r))
z['this-work-short'], dadz['this-work-short'] = data[:,0], data[:,1]

data = loadtxt("figs/walls/walls_daWB-sokolowski-%04.2f-%05.3f.dat" % (eta, delta_r))
z['sokolowski'], dadz['sokolowski'] = data[:,0], data[:,1]

data = loadtxt("figs/walls/walls_daWB-fischer-%04.2f-%05.3f.dat" % (eta, delta_r))
z['fischer'], dadz['fischer'] = data[:,0], data[:,1]

plt.figure(figsize=(6,4))

for version in versions:
  plt.plot(z[version], dadz[version], styles.plot[version], label=styles.title[version])


plt.xlim([-.5,6])
plt.legend(loc='upper right').draw_frame(False)

plt.xlabel('$z/R$')
plt.ylabel(r'$dF_1/dz$')
plt.yticks([])

if 2 < delta_r and delta_r < 2.01:
  plt.title('$\Phi(r) = \delta(\sigma - r + \delta)$     $\eta = %g$' % eta)
else:
  plt.title('$\Phi(r) = \delta(%g\sigma - r)$     $\eta = %g$' % (delta_r/2, eta))

plt.tight_layout()

savefig("figs/dadz-%d-%d.pdf" % (int(eta*10), int(delta_r)))

plt.show()
