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

if 'tmi' in sys.argv:
    moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g-tmi-dos' % (ww,ff,lenx,lenyz)
else:
    moviedir = 'figs/movies/lv/ww%.2f-ff%.2f-%gx%g-dos' % (ww,ff,lenx,lenyz)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
minlndos = 1e100
maxlndos = -1e100
numframes = 0

dataformat = 'data/lv/ww%.2f-ff%.2f-%gx%g-movie/%%06d' % (ww,ff,lenx,lenyz)
if 'tmi' in sys.argv:
    dataformat = 'data/lv/ww%.2f-ff%.2f-%gx%g-tmi-movie/%%06d' % (ww,ff,lenx,lenyz)

for frame in xrange(100000):
    basename = dataformat % frame
    try:
        e, lndos = readandcompute.e_lndos(basename)
    except:
        break
    numframes = frame+1
    minlndos = min(minlndos, lndos.min())
    maxlndos = max(maxlndos, lndos.max())
    for i in range(len(lndos)-1):
        if lndos[i] == lndos[i+1]:
            pass
        else:
            maxe = max(maxe, e[i])
            break
    for i in range(len(lndos)-1, 1, -1):
        if lndos[i-1] == lndos[i]:
            pass
        else:
            mine = min(mine, e[i]-20)
            break
    try:
        min_important_energy = readandcompute.min_important_energy(basename)
        mine = min(mine, min_important_energy-20)
    except:
        pass


print 'mine', mine
print 'maxe', maxe
print 'minlndos', minlndos
print 'maxlndos', maxlndos
print 'numframes', numframes

for frame in xrange(numframes):
    if frame % 10 == 0:
        print 'working on frame %d/%d' % (frame, numframes)
    plt.cla()

    basename = dataformat % frame

    e, lndos = readandcompute.e_lndos(basename)
    ax.plot(e, lndos, 'k-')
    datname = basename+'-lndos.dat'
    min_T = readandcompute.minT(datname)
    try:
        ax.axvline(-readandcompute.max_entropy_state(basename), color='r', linestyle=':')
        min_important_energy = readandcompute.min_important_energy(basename)
        ax.axvline(-min_important_energy, color='b', linestyle=':')
        ax.plot(e, (e+min_important_energy)/min_T + lndos[min_important_energy], 'g--')
    except:
        pass

    ax.set_xlabel(r'$E$')
    ax.set_ylim(minlndos, maxlndos)
    # ax.set_xlim(-5, -0.3)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'$\ln DOS$')
    # ax.legend(loc='best').get_frame().set_alpha(0.25)
    if 'tmi' in sys.argv:
        plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$ tmi' % (ww, ff, lenx, lenyz))
    else:
        plt.title(r'lv movie with $\lambda = %g$, $\eta = %g$, $%g\times %g$' % (ww, ff, lenx, lenyz))

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

duration = 5.0 # seconds

avconv = "avconv -y -r %g -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (numframes/duration, moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
