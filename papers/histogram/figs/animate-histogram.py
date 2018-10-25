#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os, glob, colors, argparse

matplotlib.rc('text', usetex=True)

import readnew
import re

parser = argparse.ArgumentParser(description='Animate the entropy')
parser.add_argument('subdir', metavar='s000', type=str,
                    help='the "seed" directory, typically s000')
parser.add_argument('periodic_whatever', metavar='periodic-w1.30-ff...', type=str,
                    help='the name of the job')
parser.add_argument('methods', metavar='METHOD', type=str, nargs='+',
                    help='the methods to animate')
parser.add_argument('--all-frames', action='store_true', help="plot every frame!")

args = parser.parse_args()
print(args)
subdirname = args.subdir
filename = args.periodic_whatever
suffixes = args.methods

print(sys.argv)
if subdirname[:len('data/')] == 'data/':
	subdirname = subdirname[len('data/'):]
	print('cutting redundant "data/" from first argument:', subdirname)

moviedir = 'figs/movies/%s/%s-hist' % (subdirname, filename)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
numframes = 0

dataformat = 'data/%s/%s-%%s-movie/%%06d' % (subdirname, filename)

lastframe = -1
alldone = False
for frame in range(1,100000):
    if alldone:
        break
    for suffix in suffixes:
        basename = dataformat % (suffix, frame)
        try:
            e, hist = readnew.e_and_total_init_histogram(basename)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            alldone = True
            break
        numframes = frame+1
        if len(e[hist != hist[-1]]) > 1:
            mine = min(mine, e[hist != hist[-1]].min()-5)
        if len(e[hist != hist[0]]):
            maxe = max(maxe, e[hist != hist[0]].max()+5)
    if numframes % 25 == 0 and frame != lastframe:
        print('counting %dth frame' % numframes)
        lastframe = frame

bestframe = sorted(glob.glob('data/%s/%s-%s-movie/*'
                             % (subdirname, filename, suffixes[0])))[-1]
# Now strip -transitions.dat from filename.
bestframe = re.sub('\-transitions.dat$', '', bestframe)
best_e, best_hist = readnew.e_and_total_init_histogram(bestframe)

print('best data is', bestframe)

maxhist = best_hist.max()
minhist = best_hist.min()
print('counted %d frames' % numframes)
print('mine', mine)
print('maxe', maxe)
print('minhist', minhist)
print('maxhist', maxhist)

skipby = 1
maxframes = 200
if numframes > maxframes and not args.all_frames:
    skipby = numframes // maxframes
    numframes = numframes // skipby
    print('only showing 1/%d of the frames' % skipby)
print('numframes', numframes)

for frame in range(1, numframes+1):
    if frame % 25 == 0:
        print('working on frame %d/%d' % (frame, numframes))
    plt.cla()
    ax.plot(best_e, best_hist, ':', color='0.5')

    for suffix_index in range(len(suffixes)):
        suffix = suffixes[suffix_index]
        basename = dataformat % (suffix, frame*skipby)

        try:
            e, hist = readnew.e_hist(basename)
            colors.plot(e, hist, method=suffix)
            datname = basename+'-transitions.dat'
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            print(e)
            pass

    ax.set_xlabel(r'$E$')
    ax.set_ylim(1.1*minhist, maxhist+5)
    ax.set_xlim(mine, maxe)
    ax.set_ylabel(r'$\textrm{Histogram}$')
    colors.legend(loc='lower right')

    fname = '%s/frame%06d.png' % (moviedir, frame)
    plt.savefig(fname)

duration = 10.0 # seconds

avconv = "avconv -y -r %g -i %s/frame%%06d.png -b 1000k %s/movie.mp4" % (numframes/duration, moviedir, moviedir)
os.system(avconv) # make the movie
print(avconv)
