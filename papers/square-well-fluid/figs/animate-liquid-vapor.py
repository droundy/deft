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
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]

print sys.argv

if 'tmi' in sys.argv:
    moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g-tmi' % (ww,ff,lenx,lenyz)
elif 'toe' in sys.argv:
    moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g-toe' % (ww,ff,lenx,lenyz)
else:
    moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g' % (ww,ff,lenx,lenyz)

os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
minhist = 1e100
maxhist = -1e100
numframes = 0

dataformat = 'data/lv/ww%.2f-ff%.2f-%gx%g-movie/%%06d' % (ww,ff,lenx,lenyz)
if 'tmi' in sys.argv:
    dataformat = 'data/lv/ww%.2f-ff%.2f-%gx%g-tmi-movie/%%06d' % (ww,ff,lenx,lenyz)
if 'toe' in sys.argv:
    dataformat = 'data/lv/ww%.2f-ff%.2f-%gx%g-toe-movie/%%06d' % (ww,ff,lenx,lenyz)

for frame in xrange(100000):
    basename = dataformat % frame
    try:
        e, hist = readandcompute.e_and_total_init_histogram(basename)
    except:
        break
    numframes = frame+1
    mine = min(mine, e.min() - 20)
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
    if frame % 10 == 0:
        print 'working on frame %d/%d' % (frame, numframes)
    plt.cla()

    basename = dataformat % frame
    old_min_T = min_T
    min_T = readandcompute.minT_from_transitions(basename)
    # e, diff = readandcompute.e_diffusion_estimate(basename)

    try:
        N = readandcompute.read_N(basename)
        ax.axvline(-readandcompute.max_entropy_state(basename), color='r', linestyle=':')
        ax.axvline(-readandcompute.min_important_energy(basename), color='b', linestyle=':')
        ax.axvline(-readandcompute.converged_state(basename+'-lndos.dat'), color='c', linestyle=':')
    except:
        pass

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
            theindex = numpy.nonzero(baseline_e == e[i])[0]
            if len(theindex) > 0:
                newstuff[i] -= baseline_init_hist[theindex[0]]
    else:
        newstuff = init_hist - baseline_init_hist
    ax.plot(e, newstuff, 'k-',
            label=r'%e additional iterations' % ((sum(init_hist) - sum(baseline_init_hist))/float(N)))

    ax.set_xlabel(r'$E$')
    ax.set_ylim(0, maxhist)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'histogram')
    ax.legend(loc='best').get_frame().set_alpha(0.25)
    if 'tmi' in sys.argv:
        plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$ tmi' % (ww, ff, lenx, lenyz))
    elif 'toe' in sys.argv:
        plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$ toe' % (ww, ff, lenx, lenyz))
    else:
        plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$' % (ww, ff, lenx, lenyz))

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

duration = 10.0 # seconds

avconv = "avconv -y -r %g -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (numframes/duration, moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
