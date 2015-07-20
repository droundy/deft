#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os

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

mine = 1e100
maxe = -1e100
minhist = 1e100
maxhist = -1e100
numframes = 0

for frame in xrange(100000):
    basename = 'data/mc/ww%.2f-ff%.2f-N%d-movie/%06d' % (ww,ff,N,frame)
    try:
        e, hist = readandcompute.e_and_total_init_histogram(basename)
    except:
        break
    numframes = frame+1
    mine = min(mine, e.min())
    maxe = max(maxe, e.max())
    minhist = min(minhist, hist.min())
    maxhist = max(maxhist, hist.max())

print 'mine', mine
print 'maxe', maxe
print 'minhist', minhist
print 'maxhist', maxhist
print 'numframes', numframes

min_T = None

for frame in xrange(numframes):
    print 'working on frame', frame
    plt.cla()

    basename = 'data/mc/ww%.2f-ff%.2f-N%d-movie/%06d' % (ww,ff,N,frame)
    old_min_T = min_T
    min_T = readandcompute.minT_from_transitions(basename)
    # e, diff = readandcompute.e_diffusion_estimate(basename)

    try:
        ax.axvline(-readandcompute.max_entropy_state(basename)/N, color='r', linestyle=':')
        ax.axvline(-readandcompute.min_important_energy(basename)/N, color='b', linestyle=':')
    except:
        pass

    #     T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
    #     ax.plot(u/N, T, 'k-')
    #     ax.set_ylim(0, 3)
    #     ax.axhline(minT, color='r', linestyle=':')

    #     e, hist = readandcompute.e_hist(basename)
    #     iterations = readandcompute.iterations(basename)
    #     ax.plot(e/N, 2.5*hist/hist.max(), 'k-', label=r'%e iterations' % (iterations))

    e, init_hist = readandcompute.e_and_total_init_histogram(basename)
    if min_T != old_min_T:
        print 'min_T goes from', old_min_T, 'to', min_T
        baseline_init_hist = init_hist
        baseline_e = e
    ax.plot(e, init_hist, 'b-',
            label=r'%e initialization iterations' % (sum(init_hist)/float(N)))
    # newstuff below is init_hist - baseline_init_hist, but takes into
    # account the fact that these arrays might not be the same size.
    # So we look up the element with the corresponding energy.
    if len(e) != len(baseline_e):
        newstuff = 1.0*init_hist
        for i in xrange(len(e)):
            theindex = numpy.where(baseline_e == e[i])
            if len(theindex) > 0:
                newstuff[i] -= baseline_init_hist[theindex[0]]
    else:
        newstuff = init_hist - baseline_init_hist
    ax.plot(e, newstuff, 'k-',
            label=r'%e additional iterations' % ((sum(init_hist) - sum(baseline_init_hist))/float(N)))

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(0, maxhist)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'histogram')
    ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'Histogram movie with $\lambda = %g$, $\eta = %g$, $N=%d$' % (ww, ff, N))

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

os.system("avconv -y -r 10 -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (moviedir, moviedir)) # make the movie
print("avconv -y -r 10 -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (moviedir, moviedir)) # make the movie
