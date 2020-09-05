#!/usr/bin/python3

#This program creates plots of Free Energy difference vs density 
#for many temperatures from data in kT*best.dat (or kT*best_tensor.dat)
#files which are generated as output data files by figs/new-melting.cpp

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-all-isotherm.py directory(with data files) --tensor(optional)

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


for kT in [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]:

  n = []
  gw = []
  fe_difference = []
  
  if args.tensor :
    #files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best_tensor.dat' % kT))) 
    files = sorted(list(glob.glob('newdata_tensor/phase-diagram/kT%.3f_n*_best_tensor.dat' % kT)))   #remove the "2" at the end of phase-diagram when done comparing new data
    #files = sorted(list(glob.glob('%s/kT%.3f_n*_best.dat' % (args.directory,  kT))))
  else :
    #files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
    files = sorted(list(glob.glob('newdata/phase-diagram/kT%.3f_n*_best.dat' % kT)))    #remove the "2" at the end of phase-diagram when done comparing new data
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
