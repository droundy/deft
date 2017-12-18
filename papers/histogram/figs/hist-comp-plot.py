#!/usr/bin/python2
from __future__ import division
import numpy as np
import matplotlib, sys
#if 'show' not in sys.argv:
#    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os, colors
from glob import glob

matplotlib.rc('text', usetex=True)

sys.path.insert(0, '../square-well-fluid/figs/')
import readandcompute

filebase = sys.argv[1]
#ww1.30-ff0.22-60x8
methods = [ '-tmmc', '-tmi', '-tmi2', '-tmi3', '-toe', '-toe2', '-toe3',
            '-vanilla_wang_landau', '-samc', '-satmmc', '-sad', '-sad3']
# For WLTMMC compatibility with LVMC
lvextra = glob('data/lv/%s-wltmmc*' % filebase)
split1 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra]
split2 = [i.split('-m', 1)[0] for i in split1]

for j in range(len(split2)):
    methods.append('-%s' %split2[j])
print 'methods are', methods

mine = 1e100
maxe = -1e100
maxhist = 0




for method in methods:
    try:
        print 'trying method', method
        basename = 'data/lv/%s%s-movie/000000' % (filebase, method)
        e, hist = readandcompute.e_and_total_init_histogram(basename)
        datname = basename + '-transitions.dat'
        min_T = readandcompute.minT(datname)

        hist_max = np.amax(hist)
        hist_norm = hist/hist_max

        plt.figure('error-at-energy-iterations')
        plt.subplot(2,1,1)
        colors.plot(e, hist_norm, method=method[1:])
        colors.legend()
    except:
        continue

plt.show()




