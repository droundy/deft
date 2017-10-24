#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print "Usage: directory"
    exit(1)

data_directory=sys.argv[1]
data_file=data_directory+"/plot.dat"
print "Removing plot file if it exists..."
os.system("rm "+data_file)
print "Creating new plot file [fuzzy-fmt]/"+data_file 
os.system("cat "+data_directory+"/*best.dat >>"+data_directory+"/plot.dat")  

thisdata = np.loadtxt(data_file)
print thisdata

densities=thisdata[:,1]
crystal_energies_per_atom = thisdata[:,2]
homogeneous_energies_per_atom = thisdata[:,3]
energy_differences_per_atom = thisdata[:,4]
crystal_energies_per_volume = thisdata[:,5]

print 'densities', densities
print 'crystal_energies_per_atom', crystal_energies_per_atom
print 'homogeneous_energies_per_atom', homogeneous_energies_per_atom
print 'energy_differences_per_atom', energy_differences_per_atom
print 'crystal_energies_per_volume', crystal_energies_per_volume

plot1=data_directory+"/cFEvsDen_plot.png"
plot2=data_directory+"/hFEvsDen_plot.png"
plot3=data_directory+"/Pressure_plot.png"

# Plot Crystal Free Energy per sphere vs Reduced Density
plt.plot(densities, crystal_energies_per_atom)
plt.title('Crystal Free Energy per sphere vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')
plt.savefig(plot1)

plt.figure()

# Plot Crystal Free Energy per sphere vs Reduced Density
plt.plot(densities, homogeneous_energies_per_atom)
plt.title('Homogeneous Free Energy per sphere vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')
plt.savefig(plot2)

plt.figure()

f =crystal_energies_per_atom
df=np.diff(f)
print "df=", df
n =densities
print n
dn=np.diff(n)
print dn
mid_n=n[0:len(n)-1]+dn/2
print mid_n
pressure = -(mid_n*mid_n)*(df/dn) #for fixed N and T
print pressure

# Plot Pressure vs Reduced Density
plt.plot(mid_n, pressure)
plt.title('Reduced Pressure vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Reduced Pressure')
plt.savefig(plot2)

plt.show()


#REF - delete later ----------------------------------------------------------

#EXAMPLE: cmdcpplotdat="cp plot.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"  
#         os.system(cmdcpplotdat)

#print thisdata[0] #prints row 0

#num_of_densities=len(thisdata)
#print num_of_densities

# Plot Crystal Free Energy per sphere vs Reduced Density
#plt.plot(data[:,0], data[:,1])
#plt.title('Crystal Free Energy per sphere vs Reduced Density')
#plt.xlabel('Reduced Density')
#plt.ylabel('Free Energy')
#plt.savefig('plotcopy1.png')  DEL/REF
#plt.savefig(plot1)

#print range(len(thisdata))
#for i in range(len(thisdata)):
#    print "i is", i
#    densities[i] = thisdata[i,1]
#    crystal_energies_per_atom[i] = thisdata[i,2]
#    homogeneous_energies_per_atom[i] = thisdata[i,3]
#    energy_differences_per_atom[i] = thisdata[i,4]
#    crystal_energies_per_volume[i] = thisdata[i,5]

#old pressure = ((df/dn)*(n[0:num_rd]+dn/2))-(f[0:num_rd]+df/2)
#old plt.plot(n[0:num_rd]+dn/2, pressure)
