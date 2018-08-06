#!/usr/bin/env python2

#from __future__ import division, print_function
#import sys, os
#import numpy as np
#from matplotlib import pyplot as plt
#from matplotlib import animation
#
## Specify the filepath.
#dir_path = os.path.dirname(os.path.realpath(__file__))
#sys.path.insert(0, dir_path+'/results/')
#
## Import the data.
#f_name = sys.argv[1]
#data = __import__('%s' % (f_name)) #'ising1-lnw'
#
#fig = plt.figure()
#ax = plt.axes(xlim=(np.min(data.E),np.max(data.E)), ylim=(-1000,1))
#line, = plt.plot([], [], lw=2)
#
## change this to use colors rather than line at some point
## so that each method gets its own color.
#
## Plot the background of each frame
#def init():
#    line.set_data([], [])
#    return line,
#
## animation function.  This is called sequentially
#def animate(i):
#  # for each method do this!
#  line.set_data(data.E, data.lndos[i])
#  return line,
#
## blit=True means only re-draw the parts that have changed.
#anim = animation.FuncAnimation(fig, animate, init_func=init,
#                               frames=len(data.t), interval=20, blit=True)
#
## REFERENCE: http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('papers/ising/movies/ising_lndos.mp4', fps=len(data.t)/10, extra_args=['-vcodec', 'libx264'] )
#
#plt.show()

import matplotlib, numpy, sys, os, time, glob, argparse
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'papers/histogram/figs/)')
print(os.getcwd())
if 'show' not in sys.argv:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt

import colors, readnew

matplotlib.rc('text', usetex=True)

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