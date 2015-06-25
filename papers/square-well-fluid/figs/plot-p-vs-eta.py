#!/usr/bin/python2
from __future__ import division
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rc('text', usetex=True)

plt.figure()

# input: ['data/homogeneous/ww%g-kT%g.dat' % (ww, kT) for ww in [1.3] for kT in xrange(1,11)]
ww = 1.3
for kT in np.arange(1.0, 11.0, 1.0):
    fname = 'data/homogeneous/ww%g-kT%g.dat' % (ww, kT)
    data = np.loadtxt(fname)
    plt.plot(data[:,0], data[:,1], label=r'$T=%g$' % kT)

plt.legend(loc='best')
plt.xlim(0, 0.53)
plt.ylim(-.1, 3.5)
plt.xlabel(r'$\eta$')
plt.ylabel(r'$p$')
plt.savefig('figs/p-vs-eta.pdf')
    
