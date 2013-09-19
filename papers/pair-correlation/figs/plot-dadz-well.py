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

data = loadtxt("figs/mc/a1/square-well-0.3-1.790.dat")
mc_z0, mc_da_dz = data[:,0], data[:,1]

data = loadtxt("figs/walls/square-well-dadz-this-work-%04.2f-1.790.dat" % eta)
tw_z0, tw_dadz = data[:,0], data[:,1]

data = loadtxt("figs/walls/square-well-dadz-sokolowski-%04.2f-1.790.dat" % eta)
s_z0, s_dadz = data[:,0], data[:,1]

plt.plot(tw_z0, tw_dadz, styles.plot['this-work'], label=styles.title['this-work'])
plt.plot(s_z0, s_dadz, styles.plot['sokolowski'], label=styles.title['sokolowski'])

plt.plot(mc_z0, mc_da_dz, styles.plot['mc'], label=styles.title['mc'])
plt.xlim([-.5,4])
plt.legend(loc='best').draw_frame(False)

savefig("figs/dadz-square-well-%d.pdf" % (int(eta*10)))

plt.show()
