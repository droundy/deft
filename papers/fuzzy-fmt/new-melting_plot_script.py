#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./new-melting_plot_script.py [directory where data stored] [temp]

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 3:
    print "Usage: directory  temperature"
    exit(1)

data_directory=sys.argv[1]
temp=sys.argv[2]
data_file=data_directory+"/plot_kT"+temp+".dat"
print "Removing plot file if it exists..."  #ASK is there a way to tell whether it exists so can avoid error message?
os.system("rm "+data_file)

#wait = raw_input("If not, press the ENTER key to continue program...")
print

print "Creating new plot file [fuzzy-fmt]/"+data_file 

os.system("cat "+data_directory+"/kT"+temp+"*best.dat >>"+data_file)   

thisdata = np.loadtxt(data_file)
print thisdata

densities=thisdata[:,1]
print densities
crystal_energies_per_atom = thisdata[:,5]
homogeneous_energies_per_atom = thisdata[:,4]
energy_differences_per_atom = thisdata[:,6]
crystal_energies_per_volume = thisdata[:,9]
#if want vol = 4*(1-fv)/reduced_density

plot1=data_directory+"/plot1_FEvsDen_kT"+temp+".png"
plot2=data_directory+"/plot2_Pressure_kT"+temp+".png"

# Plot Free Energy/atom vs Reduced Density
plt.plot(densities, crystal_energies_per_atom, 'b', label="Crystal Free Energy/atom")
plt.plot(densities, homogeneous_energies_per_atom, 'g', label="Homogeneous Free Energy/atom")
plt.title('Free Energy/atom vs Reduced Density at Fixed kT='+temp)
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy/atom')
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


#------------------------------------------------------------------------------
#NOTE: lattice_constant will be divided by gwstep   

#Do we need these in the plot file? - ASK!
 #crystal_energiesdensities = np.zeros_like(densities)  #initializing...
 #crystal_energies_per_volume = np.zeros_like(densities)
 #energy_differences = np.zeros_like(densities)
