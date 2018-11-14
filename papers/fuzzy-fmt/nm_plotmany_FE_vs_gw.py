#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./nm_plotmany_FE_vs_gw.py [filename*.dat]


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.')
parser.add_argument('filedat', metavar='datafile', type=str, nargs='*',
                    help='file with data to plot') 
args=parser.parse_args()

#markers = {
#   0.1: '^',
#   0.5: 'o',
#   0.01: '+',
#   0.05: 'v',
#}

plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')

for f in args.filedat:
   print 'Reading file', f
   thisdata = np.loadtxt(f)

   gw = thisdata[:,0]
   mean_FE = thisdata[:,1]
   mean_FE_uncertainty = thisdata[:,2]

   plt.errorbar(gw, mean_FE, mean_FE_uncertainty, fmt='o')
   
#plt.title('FE vs gw at kT %g n %g fv %g' % (kT, n, fv))
plt.title('FE vs gw')
plt.ylabel('FE')
#plt.xlabel('gw (with dx=%g mc_error=%g mc_constant=%g mc_prefactor=%g)' % (dx, mc_error, mc_constant, mc_prefactor))
plt.xlabel('gw')

plt.legend()

plt.show()
