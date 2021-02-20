#!/usr/bin/python3

#This program creates a plot of Free Energy difference vs gw at a specified 
#temperature and density from data in kT*n*alldat.dat (or kT*n*alldat_tensor.dat) files
#which are generated as output data files by figs/new-melting.cpp

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-FE_vs_gw.py --kT [temp] --n [density]   --fv [fraction of vacancies]  directory

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
parser.add_argument('--fv', metavar='fraction of vacancies', type=float,
                    help='fraction of vacancies - REQUIRED')
parser.add_argument('directory', metavar='directory', type=str,
                    help='directory with data to plot') 

args=parser.parse_args()

kT=args.kT
n=args.n
fv=args.fv

gw = []
fe_difference = []

#files = sorted(list(glob.glob('crystallization/kT%.3f_n%.3f_*alldat_tensor.dat' % (kT, n))))  
#files = sorted(list(glob.glob('newdata_tensor/phase-diagram/kT%.3f_n%.3f_fv%.2f*alldat_tensor.dat' % (kT, n, fv)))) 
files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.2f*alldat_tensor.dat' % (args.directory, kT, n, fv)))) 
for f in files:
    data = np.loadtxt(f)
    gw.append(data[3])
    fe_difference.append(data[6])
   
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')
plt.title('Free Energy difference vs gw  for kT=%g, n=%g, fv=%g Tensor' % (kT, n, fv))
plt.ylabel('FEdiff')
plt.xlabel('gw')
plt.plot(gw, fe_difference, '.-')

plt.show()


