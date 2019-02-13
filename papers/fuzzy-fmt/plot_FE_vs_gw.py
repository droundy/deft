#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./plot_FE_vs_gw.py --kT [temp] --n [density]  [OPTIONAL: --ten]

import os, glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot of FE vs gw at a specified temp and density.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', metavar='density', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--ten', action='store_true',
                    help='--ten for use tensor weight')
                    
args=parser.parse_args()



kT=args.kT
n=args.n

gw = []
fe_difference = []

if args.ten :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n%.3f_*alldat_tara.dat' % (kT, n))))    
else :
    files = sorted(list(glob.glob('crystallization/kT%.3f_n%.3f_*alldat.dat' % (kT, n))))
for f in files:
    data = np.loadtxt(f)
    gw.append(data[3])
    fe_difference.append(data[6])

#show constants on plot:
#dx=data[0,9]
#mc_error=data[0,10]
#mc_constant=data[0,14]
#mc_prefactor=data[0,15]
#fv=data[0,2]
#kT=data[0,0]
#n=data[0,1]
   
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')

#plt.title('FE vs gw at kT %g n %g fv %g' % (kT, n, fv))
plt.title('FE vs gw')
plt.ylabel('FE')
plt.xlabel('gw')
#plt.xlabel('gw (with dx=%g mc_error=%g mc_constant=%g mc_prefactor=%g)' % (dx, mc_error, mc_constant, mc_prefactor))
plt.plot(gw, fe_difference, '.-')

plt.show()


