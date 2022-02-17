#!/usr/bin/python3

#This program creates a plot of Free Energy difference vs gw at a specified 
#temperature and density from data in kT*n*alldat.dat (or kT*n*alldat_tensor.dat) files
#which are generated as output data files by figs/new-melting.cpp

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-FE_vs_gw-analysis.py --kT [temp] --n [density]   directory

import os, glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot of FEdiff vs gw at a specified temperature, density, and fraction of vacancies.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', metavar='density', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('directory', metavar='directory', type=str,
                    help='directory with data to plot') 

args=parser.parse_args()

kT=args.kT
n=args.n

fvs=(0, 1e-1, 1e-2, 1e-3, 1e-4)
seeds=(1,2,3,4,5)

for fv in fvs:
  for seed in seeds:
    gw = []
    fe_difference = []  
    files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.6f*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed)))) 
    #print ('%s/kT%.3f_n%.3f_fv%.6*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed))
    for f in files:
      data = np.loadtxt(f)
      gw.append(data[3])
      fe_difference.append(data[6])
    if fv == 0:
     plt.plot(gw, fe_difference, '.-', color='red', label='fv=%g seed%g' % (fv, seed))
    if fv == 0.1:
     plt.plot(gw, fe_difference, '.-', color='blue', label='fv=%g seed%g' % (fv, seed))
    if fv == 0.01:
     plt.plot(gw, fe_difference, '.-', color='green', label='fv=%g seed%g' % (fv, seed))
    if fv == 0.001:
     plt.plot(gw, fe_difference, '.-', color='purple', label='fv=%g seed%g' % (fv, seed))
    if fv == 0.0001:
     plt.plot(gw, fe_difference, '.-', color='orange', label='fv=%g seed%g' % (fv, seed))
    mcerror= data[12]

plt.legend(loc='best')      
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')
plt.title('Free Energy difference vs gw  for kT=%g, n=%g, mcerror=%g' % (kT, n, mcerror))
plt.ylabel('FEdiff')
plt.xlabel('gw')


plt.show()


