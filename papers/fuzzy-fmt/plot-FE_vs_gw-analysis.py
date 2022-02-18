#!/usr/bin/python3

#This program creates a plot of Free Energy difference vs gw at a specified 
#temperature and density from data in kT*n*alldat.dat (or kT*n*alldat_tensor.dat) files
#which are generated as output data files by figs/new-melting.cpp
#The program compute-phasediagram-analysis.py is used to run new-melting 
#for various temperatures and densities.

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-FE_vs_gw-analysis.py --kT [temp] --n [density] 
#                 --fv [OPTIONAL: enter 0, 1e-1, 1e-2, 1e-3, or 1e-4]  
#                 directory [data/phase-diagram-test-mcerr3]

import os, glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot of FEdiff vs gw at a specified temperature, density, and one or more fraction of vacancies for five seeds.')

parser.add_argument('--kT', metavar='temperature', type=float, required=True,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', metavar='density', type=float, required=True,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--fv', metavar='fraction of vacancies', type=float,
                    help='fv - OPTIONAL enter 0, 1e-1, 1e-2, 1e-3, or 1e-4')
parser.add_argument('directory', metavar='directory', type=str,
                    help='directory with data to plot ie. data/phase-diagram-test-mcerr3 or data/phase-diagram-test-mcerr4') 

args=parser.parse_args()

kT=args.kT
n=args.n
fv=args.fv

if fv is not None:
  doallfvs=0
else:
  doallfvs=1

fvs=(1e-1, 1e-2, 1e-3, 1e-4, 0)
seeds=(1,2,3,4,5)

if doallfvs == 1:
 for fv in fvs:
  for seed in seeds:
    gw = []
    fe_difference = []  
    files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.6f*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed)))) 
    for f in files:
      data = np.loadtxt(f)
      gw.append(data[3])
      fe_difference.append(data[6])
      if fv == 0.1:
        plt.plot(gw, fe_difference, '.-', color='blue')
      if fv == 0.01:
        plt.plot(gw, fe_difference, '.-', color='green')
      if fv == 0.001:
        plt.plot(gw, fe_difference, '.-', color='purple')
      if fv == 0.0001:
        plt.plot(gw, fe_difference, '.-', color='orange')
      if fv == 0:
        plt.plot(gw, fe_difference, '.-', color='red')
        #plt.plot(gw, fe_difference, '.-', label='fv=%g seed%g' % (fv, seed))
      mcerror= data[12]	
 plt.plot([],[], color='blue', label='fv=0.1 or 1e-1')
 plt.plot([],[], color='green', label='fv=0.01 or 1e-2')
 plt.plot([],[], color='purple', label='fv=0.001 or 1e-3')
 plt.plot([],[], color='orange', label='fv=0.0001 or 1e-4')
 plt.plot([],[], color='red', label='fv=0')
 plt.plot([],[], color='white', label='5 seeds for each fv')
 
        	
if doallfvs == 0:	
 for seed in seeds:
    gw = []
    fe_difference = []  
    files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.6f*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed)))) 
    for f in files:
      data = np.loadtxt(f)
      gw.append(data[3])
      fe_difference.append(data[6])
      if seed == 1:
        plt.plot(gw, fe_difference, '.-', color='blue')
      if seed == 2:
        plt.plot(gw, fe_difference, '.-', color='green')
      if seed == 3:
        plt.plot(gw, fe_difference, '.-', color='purple')
      if seed == 4:
        plt.plot(gw, fe_difference, '.-', color='orange')
      if seed == 5:
        plt.plot(gw, fe_difference, '.-', color='red')
      mcerror= data[12]	
 plt.plot([],[], color='white', label='fv=%g' % (fv))
 plt.plot([],[], color='blue', label='seed1')
 plt.plot([],[], color='green', label='seed2')
 plt.plot([],[], color='purple', label='seed3')
 plt.plot([],[], color='orange', label='seed4')
 plt.plot([],[], color='red', label='seed5')
 
plt.legend(loc='best')      
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')
plt.title('Free Energy difference vs gw  for kT=%g, n=%g, mcerror=%g' % (kT, n, mcerror))
plt.ylabel('FEdiff')
plt.xlabel('gw')

if doallfvs==1:
  plt.savefig('./figs/FE-vs-gw_kT%g_n%g_mcerror%g.pdf' % (kT, n, mcerror), transparent=True)
if doallfvs==0:
  plt.savefig('./figs/FE-vs-gw_kT%g_n%g_fv_%g_mcerror%g.pdf' % (kT, n, fv, mcerror), transparent=True)

plt.show()


