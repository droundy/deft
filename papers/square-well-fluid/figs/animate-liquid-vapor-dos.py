#!/usr/bin/python2
from __future__ import division
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

moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g-dos' % (ww,ff,lenx,lenyz)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
minlndos = 1e100
maxlndos = -1e100
numframes = 0

for frame in xrange(100000):
    basename = 'data/lv/ww%.2f-ff%.2f-%gx%g-movie/%06d' % (ww,ff,lenx,lenyz,frame)
    try:
        e, lndos = readandcompute.e_lndos(basename)
    except:
        break
    numframes = frame+1
    minlndos = min(minlndos, lndos.min())
    maxlndos = max(maxlndos, lndos.max())
    try:
        e, = readandcompute.e_hist(basename)
        mine = min(mine, e.min())
        maxe = max(maxe, e.max())
    except:
        continue

print 'mine', mine
print 'maxe', maxe
print 'minlndos', minlndos
print 'maxlndos', maxlndos
print 'numframes', numframes

for frame in xrange(numframes):
    if frame % 10 == 0:
        print 'working on frame %d/%d' % (frame, numframes)
    plt.cla()

    basename = 'data/lv/ww%.2f-ff%.2f-%gx%g-movie/%06d' % (ww,ff,lenx,lenyz,frame)

    e, lndos = readandcompute.e_lndos(basename)
    min_T = readandcompute.minT_from_transitions(basename)
    try:
        N = readandcompute.read_N(basename)
        ax.axvline(-readandcompute.max_entropy_state(basename)/N, color='r', linestyle=':')
        min_important_energy = readandcompute.min_important_energy(basename)
        ax.axvline(-min_important_energy/N, color='b', linestyle=':')
        eforline = e + N
        ax.plot(eforline/N - min_important_energy/N, eforline/min_T + lndos[min_important_energy], 'g--')
    except:
        pass
    ax.plot(e/N, lndos, 'k-')

    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(minlndos, maxlndos)
    ax.set_xlim(-5, -0.3)
    # ax.set_xlim(mine/N, maxe/N)
    ax.set_ylabel(r'histogram')
    # ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$' % (ww, ff, lenx, lenyz))

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

avconv = "avconv -y -r 2 -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
