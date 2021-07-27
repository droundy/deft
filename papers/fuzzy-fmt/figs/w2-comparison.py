#!/usr/bin/python3

# NOTE: Run this program from directory deft/papers/fuzzy-fmt
# with command: ./figs/w2-comparison.py  [show (to show plots)]

from __future__ import division, print_function
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

import wca_erf

import styles

kT = 2.5
sigma = 1

rmin=0.7
rmax=1.2

alpha, Xi, diameter = wca_erf.parameters(kT)

print(Xi)
dr = 0.05*Xi
r = np.linspace(rmin, diameter, int(rmax/dr))

figscale=1.5
plt.figure(figsize=(4*figscale, 3*figscale))

plt.axvline(diameter, color='k', linestyle='-')

fp_wca = np.exp(-wca_erf.V(r)/kT)*(-wca_erf.Vprime(r)/kT)
plt.plot(r, fp_wca, '-', color=styles.color_from_kT(kT), label=r"$f'_{wca}$ for $kT/\epsilon = %g$" % kT)

r = np.linspace(rmin, rmax, int(rmax/dr))
fp_erf = np.exp(-(r-alpha)**2/Xi**2)/Xi/np.sqrt(np.pi)
plt.plot(r, fp_erf, '--', color=styles.color_from_kT(kT),label=r"$f'_{erf}$ for $kT/\epsilon = %g$" % kT)

plt.axvline(alpha, color=styles.color_from_kT(kT), linestyle=':')

#kT=.1
kT=1

r = np.linspace(rmin, diameter, int(rmax/dr))
fp_wca = np.exp(-wca_erf.V(r)/kT)*(-wca_erf.Vprime(r)/kT)
plt.plot(r, fp_wca, '-', label=r"$f'_{wca}$ for $kT/\epsilon = %g$" % kT)

alpha, Xi, diameter = wca_erf.parameters(kT)
r = np.linspace(rmin, rmax, int(rmax/dr))
fp_erf = np.exp(-(r-alpha)**2/Xi**2)/Xi/np.sqrt(np.pi)
plt.plot(r, fp_erf, '--', color=styles.color_from_kT(kT), label=r"$f'_{erf}$ for $kT/\epsilon = %g$" % kT)
#plt.plot(r, fp_erf, '--', color=styles.color_from_kT(kT), label=r'approximation $kT/\epsilon = %g$' % kT)

plt.axvline(alpha, color=styles.color_from_kT(kT), linestyle=':')
plt.ylim(0)



kT=10

r = np.linspace(rmin, diameter, int(rmax/dr))
fp_wca = np.exp(-wca_erf.V(r)/kT)*(-wca_erf.Vprime(r)/kT)
plt.plot(r, fp_wca, '-', color=styles.color_from_kT(kT), label=r"$f'_{wca}$ for $kT/\epsilon = %g$" % kT)

alpha, Xi, diameter = wca_erf.parameters(kT)
r = np.linspace(rmin, rmax, int(rmax/dr))
fp_erf = np.exp(-(r-alpha)**2/Xi**2)/Xi/np.sqrt(np.pi)
plt.plot(r, fp_erf, '--', color=styles.color_from_kT(kT), label=r"$f'_{erf}$ for $kT/\epsilon = %g$" % kT)
#plt.plot(r, fp_erf, '--', color=styles.color_from_kT(kT), label=r'approximation $kT/\epsilon = %g$' % kT)

plt.axvline(alpha, color=styles.color_from_kT(kT), linestyle=':')
plt.ylim(0)




plt.xlim(rmin, rmax)
plt.ylabel("$f'$")
plt.xlabel('$r/\sigma$')
plt.legend(loc='best')

plt.savefig('figs/w2-comparison.pdf')
plt.show()
