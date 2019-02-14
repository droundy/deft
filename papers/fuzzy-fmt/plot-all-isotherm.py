#!/usr/bin/python2

#This program creates plots of Free Energy difference vs density 
#for many temperatures from data in kT*best.dat (or kT*best_tensor.dat)
#files which are generated as output data files by figs/new-melting.cpp

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-all-isotherm.py --tensor(optional)

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser(description='Plots FEdiff vs density for many temperatures')

parser.add_argument('--tensor', action='store_true',
                    help='--tensor for use tensor weight')

args=parser.parse_args()

#plt.figure()

#for kT in [.01, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]:
for kT in [0.5, 1]:
  n = []
  gw = []
  fe_difference = []
  
  if args.tensor :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best_tensor.dat' % kT)))    
  else :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
  for f in files:
      data = np.loadtxt(f)
      n.append(data[1])
      gw.append(data[3])
      fe_difference.append(data[6])

  plt.axhspan(0.0, -0.2, color='lightblue', alpha=0.15, lw=0)
  plt.title("Free Energy difference vs density")
  plt.plot(n, fe_difference, '.-', label='kT=%g' % (kT))
  plt.axhline(0, color='k')
  plt.xlabel('n')
  plt.ylabel('FEdiff')
  plt.legend()

plt.show()
