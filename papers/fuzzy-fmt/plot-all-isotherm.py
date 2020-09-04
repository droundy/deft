#!/usr/bin/python3

#This program creates plots of Free Energy difference vs density 
#for many temperatures from data in kT*best.dat (or kT*best_tensor.dat)
#files which are generated as output data files by figs/new-melting.cpp

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand python3 ./plot-all-isotherm.py directory(with data files) --tensor(optional)

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser(description='Plots FEdiff vs density for many temperatures')

parser.add_argument('--tensor', action='store_true',
                    help='--tensor for use tensor weight')
                    
#parser.add_argument('directory', metavar='directory', type=str,
#                    help='directory with data to plot') 

args=parser.parse_args()


#plt.figure()

#for kT in [.01, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]:
#for kT in [.01, .1, .2, .3, .4, .5, .6, .7, .8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]:
#for kT in [.4, .5, .6, .7, .8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2]:   #no data appears for lower than .4
#for kT in [.4, .5, .6, .7, .8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 4,6,8,10,12,14,16,18,20,40, 60, 80, 100, 120, 140, 160, 180, 200]:   #no data appears for lower than .4
#for kT in [.1, .2, .3, .4, .5, .6, .7, .8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 4,6,8,10,12,14,16,18,20,40]:  
#for kT in [.1, .2, .3, .4, .5, .6, .7, .8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 4,6,8,10,12,14,16,18,20,40, 60, 80, 100, 120, 140, 160, 180, 200]:  
for kT in [0.5, 1, 2, 5, 10, 20, 30, 40]:
#for kT in [1.9, 2]:
#for kT in [0.5, 1]:
  n = []
  gw = []
  fe_difference = []
  
  if args.tensor :
    #files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best_tensor.dat' % kT))) 
    files = sorted(list(glob.glob('newdata_tensor/phase-diagram5/kT%.3f_n*_best_tensor.dat' % kT)))   #remove the "2" at the end of phase-diagram when done comparing new data
    #files = sorted(list(glob.glob('%s/kT%.3f_n*_best.dat' % (args.directory,  kT))))
  else :
    #files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
    files = sorted(list(glob.glob('newdata/phase-diagram5/kT%.3f_n*_best.dat' % kT)))    #remove the "2" at the end of phase-diagram when done comparing new data
    #files = sorted(list(glob.glob('%s/kT%.3f_n*_best.dat' % (args.directory,  kT))))
  for f in files:
      data = np.loadtxt(f)
      n.append(data[1])
      gw.append(data[3])
      fe_difference.append(data[6])

  plt.axhspan(0.0, -25, color='lightblue', alpha=0.15, lw=0)
  plt.title("Free Energy difference vs density")
  plt.plot(n, fe_difference, '.-', label='kT=%g' % (kT))
  plt.axhline(0, color='k')
  plt.xlabel('n')
  plt.ylabel('FEdiff')
  plt.legend()

plt.ylim(-1,1)

plt.show()
