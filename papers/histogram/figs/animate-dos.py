#!/usr/bin/python2

import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os, glob, colors, argparse

matplotlib.rc('text', usetex=True)

import readnew

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

moviedir = 'figs/movies/%s/%s-dos' % (subdirname, filename)
os.system('rm -rf ' + moviedir)
assert not os.system('mkdir -p ' + moviedir)

fig, ax = plt.subplots()

mine = 1e100
maxe = -1e100
numframes = 0

dataformat = 'data/%s/%s-%%s-movie/%%06d' % (subdirname, filename)

lastframe = -1
alldone = False
for frame in range(100000):
    if alldone:
        break
    for suffix in suffixes:
        basename = dataformat % (suffix, frame)
        try:
            e, lndos = readnew.e_lndos(basename)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            alldone = True
            break
        numframes = frame+1
        mine = min(mine, e[lndos != lndos[-1]].min() - 5)
        maxe = max(maxe, e[lndos != lndos[0]].max()+5)
    if numframes % 25 == 0 and frame != lastframe:
        print('counting %dth frame' % numframes)
        lastframe = frame

bestframe = sorted(glob.glob('data/%s/%s-%s-movie/*-lndos.dat'
                             % (subdirname, filename, suffixes[0])))[-1]
best_e, best_lndos = readnew.e_lndos(bestframe)
print('best data is', bestframe)

maxlndos = best_lndos.max()
minlndos = best_lndos.min()
print('counted %d frames' % numframes)
print('mine', mine)
print('maxe', maxe)
print('minlndos', minlndos)
print('maxlndos', maxlndos)

skipby = 1
maxframes = 200
if numframes > maxframes and not args.all_frames:
    skipby = numframes // maxframes
    numframes = numframes // skipby
    print('only showing 1/%d of the frames' % skipby)
print('numframes', numframes)

for frame in range(numframes):
    if frame % 25 == 0:
        print('working on frame %d/%d' % (frame, numframes))
    plt.cla()
    ax.plot(best_e, best_lndos, ':', color='0.5')

    for suffix_index in range(len(suffixes)):
        suffix = suffixes[suffix_index]
        basename = dataformat % (suffix, frame*skipby)

        try:
            e, lndos, ps, lndostm = readnew.e_lndos_ps_lndostm(basename)
            colors.plot(e, lndos, method=suffix)
            if lndostm is not None and suffix[:2] != 'sa':
                colors.plot(e, lndostm, method=suffix+'-tm')
            datname = basename+'-lndos.dat'
            min_T = readnew.minT(datname)
            too_lo, too_hi = readnew.too_low_high_energy(datname)
            ax.axvline(-readnew.max_entropy_state(datname), color='r', linestyle=':')
            min_important_energy = int(readnew.min_important_energy(datname))
            ax.axvline(-min_important_energy, color='b', linestyle=':')
            if too_lo is not None and suffix[:3] == 'sad':
              ax.axvline(-too_lo, color='b', linestyle='--')
            if too_lo is not None and suffix[:3] == 'sad':
              ax.axvline(-too_hi, color='r', linestyle='--')
            # Uncomment the following to plot a line at the
            # min_important_energy with slope determined by min_T
            # ax.plot(e, (e+min_important_energy)/min_T + lndos[min_important_energy], colors[suffix_index]+'--')
            # ax.axvline(-readnew.converged_state(datname), color=colors.color(suffix), linestyle=':')
            # Uncomment the following to plot the lnw along with the lndos
            # e, lnw = readnew.e_lnw(basename)
            # ax.plot(e, -lnw, colors[suffix_index]+':')
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            print(e)
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
