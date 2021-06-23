#!/usr/bin/python

from __future__ import division

import sys
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

import wca_erf
import styles

#kT = 1.0
kT = 2.5
sigma = 1

alpha, Xi, diameter = wca_erf.parameters(kT)

rmin = 0.9
rmax = 1.2
dr = 0.01

figscale=1.5
plt.figure(figsize=(4*figscale, 3*figscale))

r = np.arange(rmin, diameter, dr)
plt.plot(r, wca_erf.V(r), 'k-', label='$V_{wca}$')

r = np.arange(rmin, rmax, dr)
Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))

plt.plot(r, Verf, 'g--', label=r'$V_{erf}$ for $kT/\epsilon = %g$' % kT)

plt.axhline(kT, color='g', ls=':')
plt.axvline(alpha, color='g', linestyle=':')

#kT = .1
kT = 1
r = np.arange(rmin, rmax, dr)
alpha, Xi, diameter = wca_erf.parameters(kT)
Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))
plt.plot(r, Verf, 'b--', label=r'$V_{erf}$ for $kT/\epsilon = %g$' % kT)
plt.axhline(kT, color='b', ls=':')
plt.axvline(alpha, color='b', linestyle=':')

plt.axvline(diameter, color='k', linestyle='-')

plt.xlim(rmin, rmax)
plt.ylim(0, 5)
plt.ylabel('$V/\epsilon$')
plt.xlabel('$r/\sigma$')
plt.legend(loc='best')

plt.savefig('figs/potential-plot.pdf')
plt.show()
