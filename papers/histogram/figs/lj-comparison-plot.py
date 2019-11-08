from __future__ import division, print_function
import sys, os, matplotlib
import numpy as np
from collections import OrderedDict

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 4)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

datadir = 'data/lj31/'

files = {
    'sad': {
            0.1: 'lj-sad-31-bin01',
            0.01: 'lj-sad-31-bin001',
            0.02: 'lj-sad-31-bin002',
            0.001: 'lj-sad-31-bin0001',
    },
    'wl': {
            0.1: 'lj-wl-31-bin01',
            0.01: 'lj-wl-31-bin001',
            0.02: 'lj-wl-31-bin002',
            0.005: 'lj-wl-31-bin0005',
            0.001: 'lj-wl-31-bin0001',
    },
    'inv-t-wl': {
            0.1: 'lj-inv-t-wl-31-bin01',
            0.01: 'lj-inv-t-wl-31-bin001',
            0.02: 'lj-inv-t-wl-31-bin002',
            0.005: 'lj-inv-t-wl-31-bin0005',
            0.001: 'lj-inv-t-wl-31-bin0001',
    },
    # 'samc-1e5': {
    #         0.1: 'lj-samc-31-1e5-bin01',
    #         0.01: 'lj-samc-31-1e5-bin001',
    #         0.02: 'lj-samc-31-1e5-bin002',
    #         0.005: 'lj-samc-31-1e5-bin0005',
    #         0.001: 'lj-samc-31-1e5-bin0001',
    # },
    # 'samc-1e6': {
    #         0.1: 'lj-samc-31-1e6-bin01',
    #         0.01: 'lj-samc-31-1e6-bin001',
    #         0.02: 'lj-samc-31-1e6-bin002',
    #         0.005: 'lj-samc-31-1e6-bin0005',
    #         0.001: 'lj-samc-31-1e6-bin0001',
    # },
    # 'samc-1e7': {
    #         0.1: 'lj-samc-31-1e7-bin01',
    #         0.01: 'lj-samc-31-1e7-bin001',
    #         0.02: 'lj-samc-31-1e7-bin002',
    #         0.005: 'lj-samc-31-1e7-bin0005',
    #         0.001: 'lj-samc-31-1e7-bin0001',
    # },
}

extra_files = [
        'lj-inv-t-wl-31-bin001-58',
        'lj-inv-t-wl-31-bin001-54',
        'lj-inv-t-wl-31-bin001-52',
]

plt.figure('cv-error')

for method in files:
        errors = {}
        iters = {}
        most_iters = 0
        for dE in files[method]:
                fname = files[method][dE]
                data = np.loadtxt(datadir+fname+'-cv-error.txt')
                if len(data) == 0:
                        continue
                errors[dE] = data[:,2]
                iters[dE] = data[:,0]
                if len(data[:,0]) > most_iters:
                        most_iters = len(data[:,0])
                        longest_iters = data[:,0]
        best = np.zeros(most_iters) + 1e30
        worst = np.zeros(most_iters)
        for i in range(most_iters):
                for dE in errors:
                        if i < len(errors[dE]):
                                best[i] = min(best[i], errors[dE][i])
                                worst[i] = max(worst[i], errors[dE][i])
        plt.fill_between(longest_iters, worst, best,
                         edgecolor='none', linewidth=0,
                         color=colors.color(method),
                         alpha=0.1, zorder=-51)
        colors.loglog(iters[0.01], errors[0.01], method=method)

for fname in extra_files:
        data = np.loadtxt(datadir+fname+'-cv-error.txt')
        if len(data) == 0:
                continue
        colors.loglog(data[:,0], data[:,2], method=fname)

moves = np.array([1e6, 1e13])
for i in np.arange(-8, 19, 1.0):
    colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')

# for fname in fnames:
#         method = fname # FIXME
#         data = np.loadtxt(datadir+fname+'-cv-error.txt')
#         colors.loglog(data[:,0], data[:,2], method=method)

plt.xlabel(r'$\textrm{Moves}$')
plt.ylabel(r'$\textrm{Maximum Error in }C_V/k_B$')

plt.xlim(1e7, 2e12)
plt.ylim(3e-1,1e4)

colors.legend()
plt.tight_layout()

plt.savefig('figs/' + 'lj-cv-error.pdf')

plt.figure('cv')

data = np.loadtxt(datadir+'bench-cv.txt')
colors.plot(data[:,0], data[:,1], method='bench')

for method in files:
        fname = files[method][0.01]
        data = np.loadtxt(datadir+fname+'-cv.txt')
        if len(data) == 0:
                continue
        colors.plot(data[:,0], data[:,1], method=method)

plt.xlim(data[0,0], data[-1,0])
plt.xlabel(r'$k_BT/\epsilon$')
plt.ylabel(r'$C_V/k_B$')

colors.legend()
plt.tight_layout()

plt.savefig('figs/' + 'lj-cv.pdf')

plt.show()
