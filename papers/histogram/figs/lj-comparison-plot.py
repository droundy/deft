from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np
from collections import OrderedDict

# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

datadir = 'data/lj31/'

filename_filter = None
if len(sys.argv) > 1:
        filename_filter = sys.argv[1]

fnames = [
    'tiny-lj-benchmark-0.001',

    'lj-wl-31-bin001',
    'lj-wl-31-bin002',
    'lj-wl-31-bin0005',

    'lj-inv-t-wl-31-bin001',
    'lj-inv-t-wl-31-bin002',
    'lj-inv-t-wl-31-bin0005',

    'lj-samc-31-1e5-bin001',
    'lj-samc-31-1e5-bin002',
    'lj-samc-31-1e5-bin0005',

    'lj-samc-31-1e6-bin001',
    'lj-samc-31-1e6-bin002',
    'lj-samc-31-1e6-bin0005',

    'lj-samc-31-1e7-bin001',
    'lj-samc-31-1e7-bin002',
    'lj-samc-31-1e7-bin0005',

    'lj-sad-31-bin001',
    'lj-sad-31-bin002',
]

plt.figure('cv-error')
for fname in fnames:
        method = fname # FIXME
        data = np.loadtxt(datadir+fname+'-cv-error.txt')
        colors.loglog(data[:,0], data[:,1], method=method)
plt.title('foo')
plt.xlabel('# moves')
plt.ylabel('error')
colors.legend()
if filename_filter is None:
        plt.savefig('figs/lj-cv-error.pdf')
else:
        plt.savefig('figs/lj-cv-error-%s.pdf' % (filename_filter))

plt.show()
