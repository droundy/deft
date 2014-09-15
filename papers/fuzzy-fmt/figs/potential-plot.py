#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

import styles

kT = 1.0

sigma = ( 2.0/(1 + 6*np.sqrt(np.log(2)*kT)) )**(1.0/6)
a = sigma/np.sqrt(np.pi)/(6*np.log(2) + np.sqrt(np.log(2)/kT))

rmin = 0.9*sigma

r = np.arange(rmin, 2**(1.0/6), 0.01)
plt.plot(r, (4/r**12 - 4/r**6 + 1)/36, label='WCA')

r = np.arange(rmin, 1.1*sigma, 0.01)
Verf = -kT*np.log(0.5*(erf ((r-sigma)/a)+1))
plt.plot(r, Verf, '-', label='erf')

plt.axvline(sigma)

plt.xlim(rmin, 1.1*sigma)
plt.ylabel('$V$')
plt.xlabel('$r$')
plt.legend(loc='best')

plt.savefig('figs/potential-plot.pdf')
