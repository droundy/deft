#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./nm_hist.py

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os

seeds=50

for gw in [.1]: 
#for gw in [0.01, 0.05, 0.1, 0.5, 0.7]: 
  print "gw=%g" % (gw)
  #for mcconstant in [300, 800]: 
  for mcconstant in [5, 20, 50, 100]: 
  #for mcconstant in [5, 20, 50, 100, 200, 500, 1000]: 
    print "mcconstant=%g" % (mcconstant)
    for seed in range(1, seeds+1):   
      print "seed=%g" % (seed)
      os.system('figs/new-melting.mkdat --kT 2 --n 1.3 --gw %g  --fv 0 --dx .5  --mc-error .0001   --mc-constant %g --mc-prefactor 10000 --seed %g | tail -n 1   >> Hist_%g_10000_gw%g_mcerror0.0001.dat' % (gw, mcconstant, seed, mcconstant, gw)) 
    os.system('./Histogram.py Hist_%g_10000_gw%g_mcerror0.0001.dat --gw %g --mcconstant %g --seeds %g >> Hist_10000_mcerror0.0001.dat' % (mcconstant, gw, gw, mcconstant, seeds)) 
