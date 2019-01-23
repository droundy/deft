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
hfe = np.array(hfe)
cfe = np.array(cfe)
invn = np.array(invn)
n = np.array(n)

print("length of n=%g" %(len(n)))
for i in range(len(n)-1):
  print(n[i])

# Plot Free Energy/atom vs 1/Reduced Density

functions = np.vstack((np.ones_like(invn),
                       invn**-1,
                       invn**-2,
                       invn**-3,
                       invn**-4,
                       invn**-5,
                       invn**-6)).T
pressure_functions = np.vstack((np.zeros_like(invn),
                                invn**-2,
                                2*invn**-3,
                                3*invn**-4,
                                4*invn**-5,
                                5*invn**-6,
                                6*invn**-7)).T
A = np.linalg.lstsq(functions, cfe)
coeff = A[0]
print('residuals', A[1])
print('coeff', coeff)
fit_cfe = np.dot(functions, coeff)
plt.plot(invn, fit_cfe, label="fit crystal free energy")

plt.plot(invn, hfe, 'red', label="Homogeneous Free Energy/atom")
plt.plot(invn, cfe, 'blue', label="Crystal Free Energy/atom")
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
fit_p = np.dot(pressure_functions, coeff)
plt.plot(invn, fit_p, label="fit crystal pressure")
plt.plot(mid_invn, hpressure, label="homogeneous pressure", color='r')
plt.plot(mid_invn, cpressure, label="crystal pressure", color='blue')
plt.title("Reduced Pressure vs 1/Reduced Density at Fixed kT=%g" % (kT))
plt.xlabel('1/Reduced Density')
plt.ylabel('Reduced Pressure')
plt.legend()
#plt.savefig(plot2)

plt.figure()

mid_hfe = 0.5*(hfe[1:] + hfe[:-1])
mid_cfe = 0.5*(cfe[1:] + cfe[:-1])

mid_h_gibbs = mid_hfe + mid_invn*hpressure
mid_c_gibbs = mid_cfe + mid_invn*cpressure
fit_c_gibbs = fit_cfe + invn*fit_p

# Plot Gibbs Free Energy/atom vs 1/Reduced Density
zoom_volume = 0.99
plt.plot(fit_p, fit_c_gibbs - fit_p*zoom_volume, 'b:', label="fit crystal")

# Solve for crossing point here:

# for i in range(len(1,mid_c_gibbs)):
#     for j in range(1,len(mid_h_gibbs)):
#         if hpressure[j] > cpressure[i-1]:


plt.plot(hpressure, mid_h_gibbs - hpressure*zoom_volume, 'r.-', label="Homogeneous Free Energy/atom")
plt.plot(cpressure, mid_c_gibbs - cpressure*zoom_volume, 'b.-', label="Crystal Free Energy/atom")
plt.title("Free Energy/atom vs 1/Reduced Density at Fixed kT=%g" % (kT))
plt.xlabel('$p$')
plt.ylabel('Gibbs Free Energy/atom - pressure*(%g volume unit)' % zoom_volume)
plt.legend()
#plt.savefig(plot1)


plt.show()


