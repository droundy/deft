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
plt.plot(n_homogeneous, kT, label='liquid', color='red')
plt.plot(n_crystal, kT, label='solid', color='blue')

#for i in range(len(kT)-1):
#   plt.plot([n_homogeneous[i], n_crystal[i]],[kT[i],kT[i]], color='gray', lw=2)
#   plt.plot([0, n_homogeneous[i]], [kT[i],kT[i]], color='red')
#   plt.plot([0, n_homogeneous[5]], [kT[5],kT[5]], color='red')  #FIX!
#   plt.plot([n_crystal[i], n_crystal[len(kT)-1]], [kT[i],kT[i]], color='blue')
   
plt.fill_betweenx(kT, 0, n_homogeneous, color='red')  
plt.fill_betweenx(kT, n_homogeneous, n_crystal, color='gray') 
plt.fill_betweenx(kT, n_crystal, 1, color='blue')  
   
plt.title("Temperature vs Density")
plt.legend(loc='best')
plt.xlabel('Density')
plt.ylabel('kT')

plt.figure()

#Pressure vs Temperature Phase Diagram
plt.fill_between(kT, 0*P, P, color='red')
plt.fill_between(kT, P, P+100, color='blue')
plt.plot(kT, P, color='black')
# for i in range(len(kT)):
#    plt.plot([0, kT[i]], [P[i],P[i]], color='darkblue')
#    plt.plot([kT[i], kT[len(kT)-1]], [P[i],P[i]], color='deepskyblue')
plt.ylim(0, 40)
plt.xlim(kT.min(), kT.max())
plt.title("Pressure vs Temperature")
plt.xlabel('kT')
plt.ylabel('Pressure')
plt.show()
