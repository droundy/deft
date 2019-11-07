from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np
from collections import OrderedDict

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
        sys.argv.remove('noshow')
import matplotlib.pyplot as plt
from glob import glob
import colors

datadir = 'data/lj31/'

filename_filter = None
if len(sys.argv) > 1:
        filename_filter = sys.argv[1]

fnames = [
    # 'tiny-lj-benchmark-0.001',

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
fnames = [f for f in fnames if '001' in f]

plt.figure('cv-error')

moves = np.array([1e6, 1e13])
for i in np.arange(-8, 19, 1.0):
    colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')
    print('diagonal', moves, 10**i/np.sqrt(0.1*moves))

for fname in fnames:
        method = fname # FIXME
        data = np.loadtxt(datadir+fname+'-cv-error.txt')
        colors.loglog(data[:,0], data[:,2], method=method)

plt.title(r'$\Delta E = 0.001\epsilon$')

plt.xlabel(r'$\textrm{Moves}$')
plt.ylabel(r'$\textrm{Maximum Error in }C_V$')

plt.xlim(1e7, 2e12)
plt.ylim(3e-1,4e3)

colors.legend()
plt.tight_layout()


if filename_filter is None:
        plt.savefig('figs/' + 'lj-cv-error.pdf')
else:
        plt.savefig('figs/lj-cv-error-%s.pdf' % (filename_filter))

plt.show()
