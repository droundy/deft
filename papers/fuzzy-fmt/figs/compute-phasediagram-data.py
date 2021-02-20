#!/usr/bin/python3

#This program takes data from *best_tensor.dat data files generated by figs/new-melting.cpp
#and located at deft/papers/fuzzy-fmt/newdata_tensor/phase-diagram4  (edit later )
#and creates 'figs/crystal-data/kT-%.3f.dat' data files with the 3-column format: density, hfe, cfe.

#This program is used in conjuction with ./figs/plot-phasediagram-data.py which produces
#temperature vs density, and pressure vs temperature phase diagrams as well as plots
#of P-vs-T at fixed n, P-vs_V at fixed T, P-vs-n at fixed T, and T-vs-n at fixed P
#from data stored in the 'figs/crystal-data/kT-%.3f.dat' files.

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand: ./figs/compute-phasediagram-data.py

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse
import sys

os.system('mkdir -p figs/crystal-data')

kT_data = []
density_data = []   #index corresponds to kT
pressure_data = []  #index corresponds to kT

#for kT in np.arange(0.1, 1.15, 0.05):   #data files with these temperatures will be plotted
#for kT in np.arange(0.1, 2.05, 0.05):  #original
#for kT in np.arange(0.4, 2.05, 0.05):   # new normal
#for kT in (1, 2, 4, 6, 8, 10, 12, 14, 16, 18):
for kT in (0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3): #for paper
#for kT in np.arange(0.1, 1.05, 0.05):   #data files with these temperatures will be plotted  DEBUG
										#values above and below this range do not currrently work   DEBUG
   
   n = []
   invn = []
   hfe = []
   cfe = []

   files = sorted(list(glob.glob('newdata_tensor/phase-diagram4/kT%.3f_n*_best_tensor.dat' % kT)))

   if len(files) == 0:
	   continue
   for f in files:
      data = np.loadtxt(f)
      n.append(data[1])     #density
      invn.append(1/data[1])
      hfe.append(data[4])   #homogeneous free energy/atom
      cfe.append(data[5])   #crystal free energy/atom
   hfe = np.array(hfe)
   cfe = np.array(cfe)
   n = np.array(n)
   np.savetxt('figs/crystal-data/kT-%.3f.dat' % kT, np.stack((n, hfe, cfe), axis=1), fmt='%.6g')
