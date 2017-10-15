#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt

#os.system('rm newmeltdataout.dat')
#os.system('rm plot.dat')

#firstlongrun for rdensity in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.9]:
#run2 for rdensity in [1.4, 1.42, 1.44, 1.46, 1.48, 1.5, 1.52, 1.54, 1.56, 1.58, 1.6, 1.62, 1.64, 1.66, 1.68, 1.7, 1.72, 1.74, 1.76, 1.78]:
#for rdensity in [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]:
#for rdensity in [1.4, 1.5] :
#    os.system('figs/new-melting.mkdat 2 %g -1 -1' % (rdensity)) 
#    os.system('cat newmeltbestdata.dat >> plot.dat')    

#os.system('gnuplot newmeltbestdata_plot.gnu')


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

