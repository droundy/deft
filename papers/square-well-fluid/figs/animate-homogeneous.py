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
ff = float(sys.argv[2])
#arg ff = [0.1,0.2,0.3,0.4,0.5]
N = int(sys.argv[3])
#arg N = [500]

moviedir = 'figs/movies/ww%.2f-ff%.2f-N%d' % (ww, ff, N)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

sleeptime = 30*60 # 30 minutes

for frame in xrange(100000):
    plt.cla()

    basename = 'data/mc/ww%.2f-ff%.2f-N%d' % (ww,ff,N)
    e, diff = readandcompute.e_diffusion_estimate(basename)

    try:
        ax.axvline(-readandcompute.max_entropy_state(basename)/N, color='r', linestyle=':')
        ax.axvline(-readandcompute.min_important_energy(basename)/N, color='b', linestyle=':')

        T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
        ax.plot(u/N, T, 'k-')
        ax.set_ylim(0, 3)
        ax.axhline(minT, color='r', linestyle=':')

        e, hist = readandcompute.e_hist(basename)
        iterations = readandcompute.iterations(basename)
        ax.plot(e/N, 2.5*hist/hist.max(), 'k-', label=r'%e iterations' % (iterations))
    except:
        pass

    e, init_hist = readandcompute.e_and_total_init_histogram(basename)
    ax.plot(e, 2.5*init_hist/init_hist.max(), 'b-',
            label=r'%e initialization iterations' % (sum(init_hist)))

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(ymin=0)
    ax.set_ylabel(r'$T$')
    ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'Histogram movie with $\lambda = %g$, $\eta = %g$, $N=%d$' % (ww, ff, N))

    fname = '%s/frame%03d.png' % (moviedir, frame)
    plt.savefig(fname)
    os.system("convert  -delay 50 %s/frame*.png %s/movie.gif" % (moviedir, moviedir)) # make the movie

    time.sleep(sleeptime)

plt.show()
