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
	print("usage: python %s s000 periodic-whatever toe-golden [toe tmi]?" % (sys.argv[0]))
subdirname = sys.argv[1]
filename = sys.argv[2]
suffixes = sys.argv[3:]

print sys.argv
if subdirname[:len('data/')] == 'data/':
	subdirname = subdirname[len('data/'):]
	print 'cutting redundant "data/" from first argument:', subdirname

moviedir = 'figs/movies/%s/%s-dos' % (subdirname, filename)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
minlndos = 1e100
maxlndos = -1e100
numframes = 0

dataformat = 'data/%s/%s-%%s-movie/%%06d' % (subdirname, filename)
colors = ['k', 'b', 'r', 'g', 'c', 'm']

lastframe = -1
alldone = False
for frame in xrange(100000):
    if alldone:
        break
    for suffix in suffixes:
        basename = dataformat % (suffix, frame)
        try:
            e, lndos = readandcompute.e_lndos(basename)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            alldone = True
            break
        minlndos = min(minlndos, lndos.min())
        maxlndos = max(maxlndos, lndos.max())
        numframes = frame+1
        mine = min(mine, e[lndos == lndos[-1]].max() - 20)
        maxe = max(maxe, e[lndos == lndos[0]].min())
    if numframes % 25 == 0 and frame != lastframe:
        print 'counting %dth frame' % numframes
        lastframe = frame

bestframe = ''
for frame in xrange(lastframe, 100000):
    basename = dataformat % (suffixes[0], frame)
    try:
        e, lndos = readandcompute.e_lndos(basename)
    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        break
    bestframe = basename
    if frame % 100 == 0:
        print 'examining %dth most golden' % frame
best_e, best_lndos = readandcompute.e_lndos(bestframe)
print 'best data is', bestframe

print 'counted %d frames' % numframes
print 'mine', mine
print 'maxe', maxe
print 'minlndos', minlndos
print 'maxlndos', maxlndos

skipby = 1
maxframes = 200
if numframes > maxframes:
    skipby = numframes // maxframes
    numframes = numframes // skipby
    print 'only showing 1/%d of the frames' % skipby
print 'numframes', numframes

for frame in xrange(numframes):
    if frame % 10 == 0:
        print 'working on frame %d/%d' % (frame, numframes)
    plt.cla()
    ax.plot(best_e, best_lndos, ':', color='0.5')

    for suffix_index in range(len(suffixes)):
        suffix = suffixes[suffix_index]
        basename = dataformat % (suffix, frame*skipby)

        try:
            e, lndos = readandcompute.e_lndos(basename)
            ax.plot(e, lndos, colors[suffix_index]+'-', label=suffix)
            datname = basename+'-lndos.dat'
            min_T = readandcompute.minT(datname)
            ax.axvline(-readandcompute.max_entropy_state(datname), color='r', linestyle=':')
            min_important_energy = readandcompute.min_important_energy(datname)
            ax.axvline(-min_important_energy, color='b', linestyle=':')
            ax.plot(e, (e+min_important_energy)/min_T + lndos[min_important_energy], colors[suffix_index]+'--')
            ax.axvline(-readandcompute.converged_state(datname), color=colors[suffix_index], linestyle=':')
            e, lnw = readandcompute.e_lnw(basename)
            ax.plot(e, -lnw, colors[suffix_index]+':')
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            pass

    ax.set_xlabel(r'$E$')
    ax.set_ylim(1.1*minlndos, maxlndos+5)
    # ax.set_xlim(-5, -0.3)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'$\ln DOS$')
    # ax.legend(loc='best').get_frame().set_alpha(0.25)
    plt.title(r'lv movie from %s ($T_{min} = %g$)' % (filename, min_T))
    plt.legend(loc='lower right')

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

duration = 10.0 # seconds

avconv = "avconv -y -r %g -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (numframes/duration, moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
