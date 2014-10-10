#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

import wca_erf

import styles

kT = 1
sigma = 1

alpha, Xi, diameter = wca_erf.parameters(kT)
rmin = alpha - 3*Xi
rmin = 0
rmax = alpha + 3*Xi
if rmax < diameter:
    rmax = diameter

r = np.arange(rmin, diameter, 0.01*Xi)

figscale=1.9
plt.figure(figsize=(4*figscale,3*figscale))

plt.axvline(diameter, color='k', linestyle=':')

fp_wca = np.exp(-wca_erf.V(r)/kT)*(-wca_erf.Vprime(r)/kT)
plt.plot(r, fp_wca, label="$f'$ WCA")

r = np.arange(rmin, rmax, 0.01*Xi)
fp_erf = np.exp(-(r-alpha)**2/Xi**2)/Xi/np.sqrt(np.pi)
plt.plot(r, fp_erf, '-', label=r"$f'$ $\alpha$-matching")

plt.axvline(alpha, color='k', linestyle=':')

plt.title(r"$f'$ with $kT/\epsilon=%g$" % kT)
plt.xlim(rmin, rmax)
plt.ylabel("$f'$")
plt.xlabel('$r$')
plt.legend(loc='best')

plt.savefig('figs/w2-comparison.pdf')
