#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-phasediagram.py [filename]

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a P vs T phase diagram.')

parser.add_argument('filedat', metavar='datafile', type=str,
                    help='file with data to plot') 
args=parser.parse_args()

thisdata = np.loadtxt(args.filedat)
kT=thisdata[:,0]   
n_homogeneous=thisdata[:,1]
n_crystal=thisdata[:,2]  
P=thisdata[:,3]

#print(len(kT))

#Temperature vs Density Phase Diagram
plt.plot(n_homogeneous, kT, color='darkblue')
plt.plot(n_crystal, kT, color='deepskyblue')
#for i in range(len(kT)-1):
     #plt.axvspan(0.84, 0.86, color='black', alpha=0.15, lw=1)
 #    plt.axvspan(n_homogeneous[i], n_crystal[i], color='black', alpha=0.15, lw=1)
plt.title("Temperature vs Density")
plt.xlabel('Density')
plt.ylabel('kT')

plt.figure()

#Pressure vs Temperature Phase Diagram
plt.plot(kT, P, color='black')
plt.title("Pressure vs Temperature")
plt.xlabel('kT')
plt.ylabel('Pressure')
plt.show()
