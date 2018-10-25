from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colors
import readnew

inputs = [
    #('vanilla_wang_landau',
    #   'data/s000/periodic-ww1.30-ff0.30-N500-vanilla_wang_landau-minE-movie/002000-lndos.dat'),
    ('samc-1e+06', 'data/s000/periodic-ww1.30-ff0.30-N500-samc-1e+06-movie/004000-lndos.dat'),
    #('sad3', 'data/s000/periodic-ww1.30-ff0.30-N500-sad3-movie/002000-lndos.dat'),
    #('sad3-test', 'data/s000/periodic-ww1.30-ff0.30-N500-sad3-test-movie/002000-lndos.dat')
]

for method, fname in inputs:
    print('plotting', method)
    e, lndos = readnew.e_lndos(fname)
    colors.plot(e, lndos, method=method)

emax = -34
emin = -3139
eminimportant = -2503
eSmax = -1206
Smin = -3700.72
plt.axvline(eminimportant)
plt.axvline(emax)
plt.axvline(emin)
plt.axvline(eSmax)
plt.ylabel('$S/k_B$')
plt.ylim(Smin, 0)
plt.xlim(emin, emax)
plt.xlabel('E')
colors.legend()

plt.savefig('figs/N500-lndos-comparison.pdf')
#plt.show()
