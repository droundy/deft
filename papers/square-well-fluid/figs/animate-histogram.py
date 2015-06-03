#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readandcompute

ww = float(sys.argv[1])
#arg ww = [1.3]
ffs = eval(sys.argv[2])
#arg ffs = [[0.1,0.2,0.3]]
lenx = float(sys.argv[3])
#arg lenx = [50,100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]

fig, ax = plt.subplots()

os.system('rm -f figs/movie-ww%.2f-%gx%g-frame*.png' % (ww,lenx,lenyz))

all_colors = ['b', 'r', 'g', 'k']
colors = {}
for i in xrange(len(ffs)):
    colors[ffs[i]] = all_colors[i]

sleeptime = 1 # second
maxsleeptime = 30*60 # 30 minutes

for frame in xrange(100000):
    plt.cla()
    for ff in ffs:
        basename = 'data/lv/ww%.2f-ff%.2f-%gx%g' % (ww,ff,lenx,lenyz)
        e, diff = readandcompute.e_diffusion_estimate(basename)

        N = readandcompute.read_N(basename);
        T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
        ax.plot(u/N, T, 'k-')
        ax.set_ylim(0, 3)
        ax.axhline(minT, color='r', linestyle=':')

        ax.axvline(-readandcompute.max_entropy_state(basename)/N, color=colors[ff], linestyle=':')
        ax.axvline(-readandcompute.min_important_energy(basename)/N, color=colors[ff], linestyle='--')

        e, hist = readandcompute.e_hist(basename)
        iterations = readandcompute.iterations(basename)
        ax.plot(e/N, 2.5*hist/hist.max(), colors[ff]+'-', label=r'$\nu = %g$, %7e iterations' % (ff, iterations))

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(ymin=0)
    ax.set_ylabel(r'$T$')
    ax.legend(loc='best')
    plt.title(r'Histogram and internal energy movie with $\lambda = %g$' % (ww))

    fname = 'figs/movie-ww%.2f-%gx%g-frame%03d.png' % (ww,lenx,lenyz, frame)
    print 'saving', fname
    plt.savefig(fname)
    os.system("convert  -delay 50 figs/movie-ww%.2f-%gx%g-frame*.png figs/movie-ww%.2f-%gx%g.gif"
              % (ww,lenx,lenyz, ww,lenx,lenyz)) # make the movie

    time.sleep(sleeptime)
    sleeptime *= 2
    if sleeptime > maxsleeptime:
        sleeptime = maxsleeptime

plt.show()
