from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colors
import readnew

plt.figure(figsize=(5, 4))

inputs = [
    #('vanilla_wang_landau',
    #   'data/s000/periodic-ww1.30-ff0.30-N500-vanilla_wang_landau-minE-movie/002000-lndos.dat'),
    ('samc-1e+06', 'data/s000/periodic-ww1.50-ff0.17-N256-samc-100000-movie/000600-lndos.dat'),
    #('sad3', 'data/s000/periodic-ww1.30-ff0.30-N500-sad3-movie/002000-lndos.dat'),
    #('sad3-test', 'data/s000/periodic-ww1.30-ff0.30-N500-sad3-test-movie/002000-lndos.dat')
]

for method, fname in inputs:
    print('plotting', method)
    e, lndos = readnew.e_lndos(fname)
    colors.plot(e, lndos, method=method)

emax = -34
emin = -2000
eminimportant = -915
eSmax = -509
Smin = -1500
plt.axvline(eminimportant, linestyle=':', color='tab:gray')
plt.axvline(eSmax, linestyle=':', color='tab:gray')
#plt.axvline(emin)
#plt.axvline(eSmax)
plt.ylabel('$S/k_B$')
plt.ylim(Smin, 0)
plt.xlim(emin, emax)
plt.xlabel('$E$')
colors.legend()

plt.tight_layout()
plt.savefig('figs/N256-lndos-comparison.pdf')
#plt.show()
