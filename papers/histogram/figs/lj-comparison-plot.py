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
    # 'bench': {
    #         # 0.01: 'tiny-lj-benchmark',
    #         # 0.01: 'lj-31-benchmark',
    # },
    'sad': {
            0.1: 'lj-sad-31-bin01',
            0.01: 'lj-sad-31-bin001',
            0.02: 'lj-sad-31-bin002',
            0.001: 'lj-sad-31-bin0001',
    },
    'WL $E_{\min}=-133.53\epsilon$': {
            0.1: 'lj-wl-31-bin01',
            0.01: 'lj-wl-31-bin001',
            0.02: 'lj-wl-31-bin002',
            0.005: 'lj-wl-31-bin0005',
            0.001: 'lj-wl-31-bin0001',
    },
    r'$1/t$-WL $E_{\min}=-133.53\epsilon$': {
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

extra_files = {
        r'$1/t$-WL $E_{\min}=-133.58\epsilon$': 'lj-inv-t-wl-31-bin001-58',
        r'$1/t$-WL $E_{\min}=-133.54\epsilon$': 'lj-inv-t-wl-31-bin001-54',
        r'$1/t$-WL $E_{\min}=-133.52\epsilon$': 'lj-inv-t-wl-31-bin001-52',
}

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
        if 'sad' in method:
                # plt.fill_between(longest_iters, worst, best,
                #                  edgecolor='none', linewidth=0,
                #                  color=colors.color(method),
                #                  alpha=0.1, zorder=-51)
                colors.loglog(iters[0.001], errors[0.001], method=r'SAD $\Delta E=0.001\epsilon$')
                colors.loglog(iters[0.01], errors[0.01], method=r'SAD $\Delta E=0.01\epsilon$')
                colors.loglog(iters[0.1], errors[0.1], method=r'SAD $\Delta E=0.1\epsilon$')
        else:
                colors.loglog(iters[0.01], errors[0.01], method=method)

for method in extra_files:
        fname = extra_files[method]
        data = np.loadtxt(datadir+fname+'-cv-error.txt')
        if len(data) == 0:
                continue
        colors.loglog(data[:,0], data[:,2], method=method)

moves = np.array([1e6, 1e13])
for i in np.arange(-8, 19, 1.0):
    colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')

# for fname in fnames:
#         method = fname # FIXME
#         data = np.loadtxt(datadir+fname+'-cv-error.txt')
#         colors.loglog(data[:,0], data[:,2], method=method)

plt.xlabel(r'Moves')
plt.ylabel(r'Maximum Error in $C_V/k_B$')

plt.xlim(1e7, 3e12)
plt.ylim(1e-1,1e4)

colors.legend()
plt.tight_layout()

plt.savefig('figs/' + 'lj-cv-error.pdf')

plt.figure('cv')

data = np.loadtxt(datadir+'bench-cv.txt')

wideax = plt.gca();
# inset axes....
inset_rect = [0.08, 0.43, 0.46, 0.57]
zoomax = wideax.inset_axes(inset_rect)

zoomax.set_xlim(data[0,0], data[-1,0])

zoom_ymin = 89
zoom_ymax = 125
zoomax.set_ylim(zoom_ymin, zoom_ymax)
colors.plot(data[:,0], data[:,1], method='bench', axes=zoomax)
data = np.loadtxt(datadir+'wide-cv.txt')
wideax.set_xlim(0, data[-1,0])
wide_ymin = 89
wide_ymax = 190
wideax.set_ylim(wide_ymin, wide_ymax)
colors.plot(data[:,0], data[:,1], method='bench', axes=wideax)

for method in files:
        fname = files[method][0.01]
        data = np.loadtxt(datadir+fname+'-cv.txt')
        widedata = np.loadtxt(datadir+fname+'-wide-cv.txt')
        if len(data) == 0:
                continue
        if 'sad' in method:
                colors.plot(data[:,0], data[:,1], method=r'SAD $\Delta E=0.01\epsilon$', axes=zoomax)
                colors.plot(widedata[:,0], widedata[:,1], method=r'SAD $\Delta E=0.01\epsilon$', axes=wideax)
                fname = files[method][0.001]
                if os.path.exists(datadir+fname+'-cv.txt'):
                        data = np.loadtxt(datadir+fname+'-cv.txt')
                        colors.plot(data[:,0], data[:,1], method=r'SAD $\Delta E=0.001\epsilon$', axes=zoomax)
                        widedata = np.loadtxt(datadir+fname+'-wide-cv.txt')
                        colors.plot(widedata[:,0], widedata[:,1], method=r'SAD $\Delta E=0.001\epsilon$', axes=wideax)
        else:
                colors.plot(data[:,0], data[:,1], method=method, axes=zoomax)
                colors.plot(widedata[:,0], widedata[:,1], method=method, axes=wideax)
for method in extra_files:
        fname = extra_files[method]
        if os.path.exists(datadir+fname+'-cv.txt'):
                data = np.loadtxt(datadir+fname+'-cv.txt')
                if len(data) == 0:
                        continue
                colors.plot(data[:,0], data[:,1], method=method, axes=zoomax)

plt.xlabel(r'$k_BT/\epsilon$')
plt.ylabel(r'$C_V/k_B$')

print('hello', 0 + inset_rect[0]*(widedata[-1,0]-0))
print(widedata[0,0])
wideax.plot([inset_rect[0]*widedata[-1,0],
             data[0,0],
             data[0,0]],
            [wide_ymax, zoom_ymax, 0], '-', color='k', linewidth=0.25, alpha=0.2, zorder=160)
wideax.plot([data[0,0], data[-1,0], data[-1,0], (inset_rect[0] + inset_rect[2])*widedata[-1,0]],
            [zoom_ymax, zoom_ymax, wide_ymin, wide_ymin+inset_rect[1]*(wide_ymax-wide_ymin)],
            '-', color='k', linewidth=0.25, alpha=0.2)

# axins.set_xlabel(r'$k_BT/\epsilon$')
# axins.set_ylabel(r'$C_V/k_B$')
# axins.set_xlim(0,0.4)

colors.legend(framealpha=1, loc='lower right')

plt.tight_layout()

plt.savefig('figs/' + 'lj-cv.pdf')

plt.show()
