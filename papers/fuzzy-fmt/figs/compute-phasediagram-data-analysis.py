#!/usr/bin/python3

# This program takes data from *best_tensor.dat data files generated by
# figs/new-melting.cpp and creates kT-%.3f.dat data files with the
# 3-column format: density, hfe, cfe.

# This program is used in conjuction with plot-phasediagram-data.py which 
# produces temperature vs density, and pressure vs temperature phase diagrams 
# and plots of P-vs-T at fixed n, P-vs_V at fixed T, P-vs-n at fixed T, 
# and T-vs-n at fixed P from the kT-%.3f.dat data files.

# NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#       with comand: ./figs/compute-phasediagram-data.py

from __future__ import print_function, division

import numpy as np
import os, glob
import argparse
import sys

os.system('mkdir -p figs/crystal-data') 

kT_data = []
density_data = []   #index corresponds to kT

for seed in (1, 2, 3, 4, 5):
 for kT in (0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 
  1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3): #for paper
   
   n = []
   invn = []
   hfe = []
   cfe = []
   hpressure_nm = []
   cpressure_nm = []

   #files = sorted(list(glob.glob('data/phase-diagram/kT%.3f_n*_best_tensor.dat' % kT)))
   #files = sorted(list(glob.glob('data/phase-diagram/kT%.3f_n*seed%g_best_tensor.dat' % (kT, seed))))  
   files = sorted(list(glob.glob('data/phase-diagram_fv1e-4/kT%.3f_n*seed%g_best_tensor.dat' % (kT, seed))))

   if len(files) == 0:
	   continue
   for f in files:
      data = np.loadtxt(f)
      n.append(data[1])     #density
      invn.append(1/data[1])
      hfe.append(data[4])   #homogeneous free energy/atom
      cfe.append(data[5])   #crystal free energy/atom
      hpressure_nm.append(data[17])
      cpressure_nm.append(data[18])
   hfe = np.array(hfe)
   cfe = np.array(cfe)
   hpressure_nm = np.array(hpressure_nm)
   cpressure_nm = np.array(cpressure_nm)
   n = np.array(n)
   #np.savetxt('figs/crystal-data/kT-%.3f.dat' % kT, np.stack((n, hfe, cfe, hpressure_nm, cpressure_nm), axis=1), fmt='%.6g')
   #np.savetxt('figs/crystal-data/kT-%.3f-seed%g.dat' % (kT, seed), np.stack((n, hfe, cfe, hpressure_nm, cpressure_nm), axis=1), fmt='%.6g')
   np.savetxt('figs/crystal-data/kT-%.3f-fv1e-4-seed%g.dat' % (kT, seed), np.stack((n, hfe, cfe, hpressure_nm, cpressure_nm), axis=1), fmt='%.6g')  
   #np.savetxt('figs/crystal-data/kT-%.3f-seed%g.dat' % (kT, seed), np.stack((n, hfe, cfe), axis=1), fmt='%.6g')  