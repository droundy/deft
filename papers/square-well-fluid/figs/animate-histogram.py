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

moviedir = 'figs/movies/ww%.2f-%gx%g' % (ww,lenx,lenyz)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

all_colors = ['b', 'r', 'g', 'k']
colors = {}
for i in xrange(len(ffs)):
    colors[ffs[i]] = all_colors[i]

sleeptime = 30*60 # 30 minutes

for frame in xrange(100000):
    plt.cla()
    for ff in ffs:
        basename = 'data/lv/ww%.2f-ff%.2f-%gx%g' % (ww,ff,lenx,lenyz)
        e, diff = readandcompute.e_diffusion_estimate(basename)

        N = readandcompute.read_N(basename);
        try:
            ax.axvline(-readandcompute.max_entropy_state(basename)/N, color=colors[ff], linestyle=':')
            ax.axvline(-readandcompute.min_important_energy(basename)/N, color=colors[ff], linestyle=':')

            T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
            ax.plot(u/N, T, 'k-')
            ax.set_ylim(0, 3)
            ax.axhline(minT, color='r', linestyle=':')

            e, hist = readandcompute.e_hist(basename)
            iterations = readandcompute.iterations(basename)
            ax.plot(e/N, 2.5*hist/hist.max(), colors[ff]+'-', label=r'$\eta = %g$, %e iterations' % (ff, iterations))
        except:
            pass

        e, init_hist = readandcompute.e_and_total_init_histogram(basename)
        ax.plot(e, 2.5*init_hist/init_hist.max(), colors[ff]+'--',
                label=r'$\eta = %g$, %e initialization iterations' % (ff, sum(init_hist)))

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(ymin=0)
    ax.set_ylabel(r'$T$')
    ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'Histogram and internal energy movie with $\lambda = %g$' % (ww))

    fname = '%s/frame%03d.png' % (moviedir, frame)
    plt.savefig(fname)
    os.system("convert  -delay 50 %s/frame*.png %s/movie.gif" % (moviedir, moviedir)) # make the movie
    os.system("avconv -y -r 10 -i %s/frame%03d.png -b 1000k %s/movie.mp4" % (moviedir, moviedir)) # make the movie

    time.sleep(sleeptime)

plt.show()
