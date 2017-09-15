#!/usr/bin/python2
from __future__ import division
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os

matplotlib.rc('text', usetex=True)

import readandcompute

if len(sys.argv) < 4:
    print("Usage: python {} lv ww1.30-ff0.22-100x10 tmi toe ...".format(sys.argv[0]))
    exit(1)

subdirname = sys.argv[1]
filename = sys.argv[2]
suffixes = sys.argv[3:]

print sys.argv

moviedir = 'figs/movies/%s/%s-hist' % (subdirname, filename)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
maxhist = 0
numframes = 0

dataformat = 'data/%s/%s-%%s-movie/%%06d' % (subdirname, filename)
colors = ['k', 'b', 'r', 'g']

lastframe = -1
for frame in xrange(100000):
    for suffix in suffixes:
        basename = dataformat % (suffix, frame)
        try:
            e, hist = readandcompute.e_and_total_init_histogram(basename)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            break
        numframes = frame+1
        maxhist = max(maxhist, hist.max())
        mine = min(mine, e.min() - 20)
        maxe = max(maxe, e.max())
    if numframes % 25 == 0 and frame != lastframe:
        print 'counting %dth frame' % numframes
        lastframe = frame

print 'mine', mine
print 'maxe', maxe
print 'maxhist', maxhist

skipby = 1
if numframes > 200:
    skipby = numframes // 200
    numframes = numframes // skipby
    print 'only showing 1/%d of the frames' % skipby

print 'numframes', numframes

for frame in xrange(numframes):
    if frame % 10 == 0:
        print 'working on frame %d/%d' % (frame, numframes)
    plt.cla()

    for suffix_index in range(len(suffixes)):
        suffix = suffixes[suffix_index]
        basename = dataformat % (suffix, frame*skipby)

        try:
            e, hist = readandcompute.e_and_total_init_histogram(basename)
            ax.plot(e, hist, colors[suffix_index]+'-', label=suffix)
            datname = basename+'-transitions.dat'
            min_T = readandcompute.minT(datname)
            ax.axvline(-readandcompute.max_entropy_state(basename), color='r', linestyle=':')
            min_important_energy = readandcompute.min_important_energy(basename)
            #ax.axvline(-min_important_energy, color='b', linestyle=':')
            #ax.axvline(-readandcompute.converged_state(datname), color=colors[suffix_index], linestyle=':')
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            pass

    ax.set_xlabel(r'$E$')
    ax.set_ylim(0, maxhist)
    # ax.set_xlim(-5, -0.3)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'histogram')
    # ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'lv movie from %s ($T_{min} = %g$)' % (filename, min_T))
    plt.legend(loc='best')

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

duration = 10.0 # seconds

avconv = "avconv -y -r %g -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (numframes/duration, moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
