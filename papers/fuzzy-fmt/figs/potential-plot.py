#!/usr/bin/python3

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

kT = 2.5
sigma = 1

alpha, Xi, diameter = wca_erf.parameters(kT)

rmin = 0.9
rmax = 1.2
dr = 0.01

figscale=1
plt.figure(figsize=(4*figscale, 3*figscale))

r = np.arange(rmin, diameter, dr)
plt.plot(r, wca_erf.V(r), 'k-', label='$V_{wca}$')

r = np.arange(rmin, rmax, dr)
Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))

plt.plot(r, Verf, '--', color=styles.color_from_kT(kT), linewidth=1, label=r'$V_{erf}$  $T^* = %g$' % kT)

#plt.axhline(kT, color='g', ls=':')
#plt.axvline(alpha, color='g', linestyle=':')
plt.axvline(alpha, color=styles.color_from_kT(kT), linestyle=':')

kT = 1
r = np.arange(rmin, rmax, dr)
alpha, Xi, diameter = wca_erf.parameters(kT)
Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))
plt.plot(r, Verf, '--', color=styles.color_from_kT(kT), linewidth=1, label=r'$V_{erf}$  $T^* = %g$' % kT)
# plt.axhline(kT, color='b', ls=':')
plt.axvline(alpha, color=styles.color_from_kT(kT), linestyle=':')


# kT = 10
# r = np.arange(rmin, rmax, dr)
# alpha, Xi, diameter = wca_erf.parameters(kT)
# Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))
# plt.plot(r, Verf, color=styles.color_from_kT(kT), label=r'$V_{erf}$  $T* = %g$' % kT)
# plt.axvline(alpha, color=styles.color_from_kT(kT),  linestyle=':')

plt.tight_layout()
plt.axvline(diameter, color='k', linestyle='-')
plt.xlim(rmin, rmax)
plt.xticks(fontsize=8)
plt.ylim(0, 5)
plt.yticks(fontsize=8)
plt.ylabel('$V/\epsilon$', fontsize=10)
plt.xlabel('$r/\sigma$', fontsize=10)
plt.legend(loc='best', fontsize=9)

plt.savefig('figs/potential-plot.pdf')
plt.show()
