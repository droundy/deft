#!/usr/bin/python

from __future__ import division
import matplotlib
import sys
import styles # for our style choices for these plots
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " eta")
    exit(1)

from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

eta = float(sys.argv[1])
#arg eta = [.2,.3]
able_to_read_file = True

plt.title('$da/dz,$ $\eta = %g,$  square well' %(eta))

data = loadtxt("figs/mc/a1/inverse-sixth-0.3-rmax-5.dat")
mc_z0, mc_da_dz = data[:,0], data[:,1]

data = loadtxt("figs/walls/inverse-sixth-dadz-this-work-%04.2f-rmax-5.dat" % eta)
tw_z0, tw_dadz = data[:,0], data[:,1]

datamc = loadtxt("figs/walls/inverse-sixth-dadz-this-work-mc-%04.2f-rmax-5.dat" % eta)
twmc_z0, twmc_dadz = datamc[:,0], datamc[:,1]

data = loadtxt("figs/walls/inverse-sixth-dadz-sokolowski-%04.2f-rmax-5.dat" % eta)
s_z0, s_dadz = data[:,0], data[:,1]

plt.figure(figsize=(6,4))

plt.plot(mc_z0[::10], mc_da_dz[::10], styles.plot['mc'], label=styles.title['mc'])
plt.plot(tw_z0, tw_dadz, styles.plot['this-work'], label=styles.title['this-work'])
plt.plot(twmc_z0, twmc_dadz, styles.plot['this-work-mc'], label=styles.title['this-work-mc'])
plt.plot(s_z0, s_dadz, styles.plot['sokolowski'], label=styles.title['sokolowski'])

plt.xlim([-.5, 6])
plt.legend(loc='best').draw_frame(False)

plt.xlabel('$z/R$')
plt.ylabel(r'$dF_1/dz$')
plt.yticks([])
#plt.xticks([0,1,2,3])

plt.title('$\Phi(r) = \Theta(2.5\sigma - r)/r^6$     $\eta = %g$' % eta)

plt.tight_layout()

savefig("figs/dadz-inverse-sixth-%d.pdf" % (int(eta*10)))

plt.show()
