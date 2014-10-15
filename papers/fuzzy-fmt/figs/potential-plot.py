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

kT = 1.0
sigma = 1

alpha, Xi, diameter = wca_erf.parameters(kT)

rmin = 0.9*alpha
rmin = 0.6*sigma
rmax = diameter

r = np.arange(rmin, diameter, 0.01)
plt.plot(r, wca_erf.V(r), label='WCA')

r = np.arange(rmin, diameter, 0.01)
Verf = -kT*np.log(0.5*(erf ((r-alpha)/Xi)+1))
plt.plot(r, Verf, '-', label=r'erf $\alpha$-matching')
plt.axvline(alpha)

plt.xlim(rmin, rmax)
plt.ylabel('$V$')
plt.xlabel('$r$')
plt.legend(loc='best')

plt.savefig('figs/potential-plot.pdf')
