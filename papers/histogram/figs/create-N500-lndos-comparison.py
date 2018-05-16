from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colors
import readnew

inputs = [
    ('vanilla_wang_landau',
       'data/s000/periodic-ww1.30-ff0.30-N500-vanilla_wang_landau-movie/000509-lndos.dat'),
    ('golden', 'data/s000/periodic-ww1.30-ff0.30-N500-tmmc-golden-movie/000324-lndos.dat'),
    ('sad3', 'data/s000/periodic-ww1.30-ff0.30-N500-sad3-movie/001952-lndos.dat'),
]

for method, fname in inputs:
    print('plotting', method)
    e, lndos = readnew.e_lndos(fname)
    colors.plot(e, lndos, method=method)

colors.legend()

plt.savefig('figs/N500-lndos-comparison.pdf')
#plt.show()
