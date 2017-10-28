#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./new-melting_miniplot_script.py [directory where data stored] [temp]
#to create plots from plot.dat files already in the data directory

import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches

if len(sys.argv) < 3:
    print "Usage: directory temperature"
    exit(1)

data_directory=sys.argv[1]
temp=sys.argv[2]
data_file=data_directory+"/plot_kT"+temp+".dat"

thisdata = np.loadtxt(data_file)

densities=thisdata[:,1]
crystal_energies_per_atom = thisdata[:,2]
homogeneous_energies_per_atom = thisdata[:,3]
energy_differences_per_atom = thisdata[:,4]
crystal_energies_per_volume = thisdata[:,5]

plot1=data_directory+"/plot1_FEvsDen_kT"+temp+".png"
plot2=data_directory+"/plot2_Pressure_kT"+temp+".png"

# Plot Free Energy/atom vs Reduced Density
plt.plot(densities, crystal_energies_per_atom, 'b', label="Crystal Free Energy/atom")
plt.plot(densities, homogeneous_energies_per_atom, 'g', label="Homogeneous Free Energy/atom")
plt.title('Free Energy/atom vs Reduced Density at Fixed kT='+temp)
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')
plt.legend()
plt.savefig(plot1)

plt.figure()

f =crystal_energies_per_atom
#print "f=", f
df=np.diff(f)  #Caution: depends on order of data files!
#print "df=", df
n =densities
#print "n=", n
dn=np.diff(n)  #Caution: depends on order of data files!
#print "dn=", dn
mid_n=n[0:len(n)-1]+dn/2
#print "mid_n=", mid_n
pressure = -(mid_n*mid_n)*(df/dn) #for fixed N and T
#print "pressure =", pressure

# Plot Pressure vs Reduced Density
plt.plot(mid_n, pressure, color='red')
plt.title('Reduced Pressure vs Reduced Density at Fixed kT='+temp)
plt.xlabel('Reduced Density')
plt.ylabel('Reduced Pressure')
plt.savefig(plot2)

plt.show()

