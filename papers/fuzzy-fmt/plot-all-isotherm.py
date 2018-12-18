#!/usr/bin/python2
#From directory deft/papers/fuzzy-fmt run command ./plot-all-isotherm.py
from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

#parser = argparse.ArgumentParser(description='Creates data for a FE vs gw plot.')

#parser.add_argument('--kT', metavar='temperature', type=float,
#                    help='reduced temperature - REQUIRED')

#args=parser.parse_args()

#kT=args.kT

plt.figure()

for kT in [.01, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3]:
#for kT in [.3]:
  n = []
  gw = []
  fe_difference = []

  files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
  for f in files:
      data = np.loadtxt(f)
      n.append(data[1])
      gw.append(data[3])
      fe_difference.append(data[6])

  #plt.figure()
  #plt.plot(n, gw, '.-')
  #plt.xlabel('n')
  #plt.ylabel('gw')

  #plt.figure()
  plt.axhspan(0.0, -0.2, color='lightblue', alpha=0.15, lw=0)
  plt.title("Isotherms at mc_error=0.001, mc_constant=5, mc_pref=50000, dx=0.5")
  plt.plot(n, fe_difference, '.-', label='kT=%g' % (kT))
  plt.axhline(0, color='k')
  plt.xlabel('n')
  plt.ylabel('FED')
  plt.legend()

plt.show()
