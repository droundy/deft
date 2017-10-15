#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt

os.system('rm newmeltdataout.dat')
os.system('rm plot.dat')

for rdensity in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]:
    os.system('figs/new-melting.mkdat 2 %g -1 -1' % (rdensity)) 
    os.system('cat newmeltbestdata.dat >> plot.dat')    

# Plot Crystal Free Energy per sphere vs Reduced Density
data = np.loadtxt('plot.dat')

plt.plot(data[:,0], data[:,1])

plt.title('Crystal Free Energy per sphere vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')

plt.show()

f=data[:,3]
df=np.diff(data[:,3])
n=data[:,0]
dn=np.diff(data[:,0])
pressure = ((df/dn)*(n[0:3]+dn/2))-(f[0:3]+df/2)

# Plot Crystal Free Energy per volume vs Reduced Density
plt.plot(n[0:3]+dn/2, pressure)

plt.title('Pressure vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Pressure')

plt.show()

