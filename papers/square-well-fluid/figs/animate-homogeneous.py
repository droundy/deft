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
ffs = eval(sys.argv[2])
Ns = eval(sys.argv[3])

fig, ax = plt.subplots()

os.system('rm -f figs/movie-ww%.2f-frame*.png' % (ww))

all_colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y']
colors = {}
for i in xrange(len(ffs)):
    colors[ffs[i]] = all_colors[i]

sleeptime = 1 # second
maxsleeptime = 30*60 # 30 minutes

for frame in xrange(100000):
    plt.cla()
    for ff in ffs:
        for N in Ns:
            basename = 'data/mc/ww%.2f-ff%.2f-N%d' % (ww,ff,N)
            e, diff = readandcompute.e_diffusion_estimate(basename)

            T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
            ax.plot(u/N, T, 'k-')
            ax.set_ylim(0, 3)
            ax.axhline(minT, color='r', linestyle=':')

            ax.axvline(-readandcompute.max_entropy_state(basename)/N, color=colors[ff], linestyle=':')
            ax.axvline(-readandcompute.min_important_energy(basename)/N, color=colors[ff], linestyle='--')

            e, hist = readandcompute.e_hist(basename)
            iterations = readandcompute.iterations(basename)
            ax.plot(e/N, 2.5*hist/hist.max(), colors[ff]+'-',
                    label=r'$\eta = %g$, $N=%d$, %.0e iterations' % (ff, N, iterations))

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(ymin=0)
    ax.set_ylabel(r'$T$')
    ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'Histogram and internal energy movie with $\lambda = %g$' % (ww))

    fname = 'figs/movie-ww%.2f-frame%06d.png' % (ww, frame)
    print 'saving', fname
    plt.savefig(fname)
    os.system("convert -delay 100 figs/movie-ww%.2f-frame*.png figs/movie-ww%.2f.gif"
              % (ww, ww)) # make the movie

    time.sleep(sleeptime)
    sleeptime *= 2
    if sleeptime > maxsleeptime:
        sleeptime = maxsleeptime

plt.show()
