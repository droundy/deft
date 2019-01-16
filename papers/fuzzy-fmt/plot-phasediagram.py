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

#a=0
#b=-2
#if b<0 :
#  print b
#else : print a



#Temperature vs Density Phase Diagram
plt.plot(n_homogeneous, kT, color='red')
plt.plot(n_crystal, kT, color='deepskyblue')
for i in range(len(kT)-1):
   plt.plot([n_homogeneous[i], n_crystal[i]],[kT[i],kT[i]], color='gray', lw=2)
   plt.plot([0, n_homogeneous[i]], [kT[i],kT[i]], color='red')
   plt.plot([0, n_homogeneous[5]], [kT[5],kT[5]], color='red')  #FIX!
   plt.plot([n_crystal[i], n_crystal[len(kT)-1]], [kT[i],kT[i]], color='deepskyblue')
plt.title("Temperature vs Density")
plt.xlabel('Density     (red=liquid, light blue=crystal)')
plt.ylabel('kT')

plt.figure()

#Pressure vs Temperature Phase Diagram
plt.plot(kT, P, color='black')
for i in range(len(kT)):
   plt.plot([0, kT[i]], [P[i],P[i]], color='deepskyblue')
   plt.plot([kT[i], kT[len(kT)-1]], [P[i],P[i]], color='red')
plt.title("Pressure vs Temperature")
plt.xlabel('kT     (red=liquid, light blue=crystal)')
plt.ylabel('Pressure')
plt.show()
