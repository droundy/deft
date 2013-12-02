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
eta = float(sys.argv[1])
#arg eta = [.2,.3]
able_to_read_file = True

plt.title('$da/dz,$ $\eta = %g,$  square well' %(eta))

data = loadtxt("figs/mc/a1/inverse-sixth-0.3-rmax-5.dat")
mc_z0, mc_da_dz = data[:,0], data[:,1]

data = loadtxt("figs/walls/inverse-sixth-dadz-this-work-%04.2f-rmax-5.dat" % eta)
tw_z0, tw_dadz = data[:,0], data[:,1]

data = loadtxt("figs/walls/inverse-sixth-dadz-sokolowski-%04.2f-rmax-5.dat" % eta)
s_z0, s_dadz = data[:,0], data[:,1]

plt.figure(figsize=(5,3))

plt.plot(tw_z0/2, tw_dadz, styles.plot['this-work'], label=styles.title['this-work'])
plt.plot(s_z0/2, s_dadz, styles.plot['sokolowski'], label=styles.title['sokolowski'])

plt.plot(mc_z0[::10]/2, mc_da_dz[::10], styles.plot['mc'], label=styles.title['mc'])
plt.xlim([-.5/2,6/2.])
plt.legend(loc='best').draw_frame(False)

plt.xlabel('$z/\sigma$')
plt.ylabel(r'$\frac{da_1}{dz}$')
plt.yticks([])
plt.xticks([0,1,2,3])

plt.title('$\Phi(r) = \Theta(2.5\sigma - r)/r^6$     $\eta = %g$' % eta)

plt.tight_layout()

savefig("figs/dadz-inverse-sixth-%d.pdf" % (int(eta*10)))

plt.show()
