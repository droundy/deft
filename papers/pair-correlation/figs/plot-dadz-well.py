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

data = loadtxt("figs/mc/a1/square-well-0.3-1.790.dat")
mc_z0, mc_da_dz = data[:,0], data[:,1]

# data = loadtxt("figs/walls/square-well-dadz-this-work-%04.2f-1.790.dat" % eta)
# tw_z0, tw_dadz = data[:,0], data[:,1]

# datamc = loadtxt("figs/walls/square-well-dadz-this-work-mc-%04.2f-1.790.dat" % eta)
# twmc_z0, twmc_dadz = datamc[:,0], datamc[:,1]

data = loadtxt("figs/walls/square-well-dadz-this-work-short-%04.2f-1.790.dat" % eta)
tws_z0, tws_dadz = data[:,0], data[:,1]

data = loadtxt("figs/walls/square-well-dadz-sokolowski-%04.2f-1.790.dat" % eta)
s_z0, s_dadz = data[:,0], data[:,1]

scaleme=0.8
plt.figure(figsize=(6*scaleme,4*scaleme))

plt.plot(mc_z0[::10], mc_da_dz[::10], styles.plot['mc'], label=styles.title['mc'])
# plt.plot(tw_z0, tw_dadz, styles.plot['this-work'], label=styles.title['this-work'])
# plt.plot(twmc_z0, twmc_dadz, styles.plot['this-work-mc'], label=styles.title['this-work-mc'])
plt.plot(tws_z0, tws_dadz, styles.plot['this-work-short'], label=styles.title['this-work-short'])
plt.plot(s_z0, s_dadz, styles.plot['sokolowski'], label=styles.title['sokolowski'])

plt.xlim([-.5, 6])
plt.legend(loc='best').draw_frame(False)

plt.xlabel('$z/R$')
plt.ylabel(r'$dF_1/dz$')
plt.yticks([])
# plt.xticks([0,1,2,3])

plt.title('$\Phi(r) = \Theta(1.79\sigma - r)$     $\eta = %g$' % eta)

plt.tight_layout()

savefig("figs/dadz-square-well-%d.pdf" % (int(eta*10)))

plt.show()
