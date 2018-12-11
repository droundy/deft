#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-pressure.py --kT [temp] 

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser(description='Creates data for a FE vs gw plot.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')

args=parser.parse_args()

kT=args.kT

n = []
invn = []
hfe = []
cfe = []

files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))
for f in files:
    data = np.loadtxt(f)
    n.append(data[1])     #density
    invn.append(1/data[1])
    hfe.append(data[4])   #homogeneous free energy/atom
    cfe.append(data[5])   #crystal free energy/atom

#print("length of n=%g" %(len(n)))
#for i in range(len(n)-1):    
#  print(n[i])

# Plot Free Energy/atom vs 1/Reduced Density
plt.plot(invn, hfe, 'deepskyblue', label="Homogeneous Free Energy/atom")
plt.plot(invn, cfe, 'darkblue', label="Crystal Free Energy/atom")
plt.title("Free Energy/atom vs 1/Reduced Density at Fixed kT=%g" % (kT))
plt.xlabel('1/Reduced Density')
plt.ylabel('Free Energy/atom')
plt.legend()
#plt.savefig(plot1)

plt.figure()

dhfe=np.diff(hfe)  #Caution: depends on order of data files!
dcfe=np.diff(cfe)  #Caution: depends on order of data files!
dinvn=np.diff(invn)  #Caution: depends on order of data files!
mid_invn=invn[0:len(invn)-1]+dinvn/2
hpressure = -(dhfe/dinvn) #for fixed N and Te   
cpressure = -(dcfe/dinvn) #for fixed N and Te  

# Plot Pressure vs 1/Reduced Density
plt.plot(mid_invn, hpressure, label="homogeneous pressure", color='deepskyblue')
plt.plot(mid_invn, cpressure, label="crystal pressure", color='darkblue')
plt.title("Reduced Pressure vs 1/Reduced Density at Fixed kT=%g" % (kT))
plt.xlabel('1/Reduced Density')
plt.ylabel('Reduced Pressure')
plt.legend()
#plt.savefig(plot2)

plt.show()


