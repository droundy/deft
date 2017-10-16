#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#for arg in sys.argv[1:]:
#    print arg

if len(sys.argv) < 5:
    print "Usage: T dfv dgw [DENSITIES TO COMPUTE]..."
    exit(1)

temp=sys.argv[1]    
fvstep=sys.argv[2]
gwstep=sys.argv[3]   #lattice_constant will be divided by this number
rdensities=sys.argv[4:]  #if enter by command line and not by hand
#print rdensities
temp_num=float(temp)

os.system('rm newmeltdataout.dat')
os.system('rm plot.dat')

#rdensities=[0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]  #if enter by hand and not by command line
num_rd=len(rdensities)-1 

everything_file="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"
best_file="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"

densities = np.array([float(d) for d in rdensities])
crystal_energies = np.zeros_like(densities)
crystal_energies_per_volume = np.zeros_like(densities)
energy_differences = np.zeros_like(densities)

for i in range(len(densities)):
    rdensity_num=densities[i]
    print 'running with', rdensity_num
#   os.system('figs/new-melting.mkdat 2 %g -1 -1' % (rdensity))    #if enter by hand and not by command line
    os.system('figs/new-melting.mkdat 2 %g -1 -1' % rdensity_num)  #if enter by command line and not by hand
    os.system('cat newmeltbestdata.dat >> plot.dat')
    thisdata = np.loadtxt('newmeltbestdata.dat')
    crystal_energies[i] = thisdata[1]
    energy_differences[i] = thisdata[2]
    crystal_energies_per_volume[i] = thisdata[3]

print 'densities', densities
print 'crystal_energies', crystal_energies
print 'energy_differences', energy_differences
print 'crystal_energies_per_volume', crystal_energies_per_volume

cmdcpdat="cp newmeltdataout.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_newmeltdataout.dat" 
os.system(cmdcpdat)
cmdcpplotdat="cp plot.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"  
os.system(cmdcpplotdat)

plot1="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_FEvsDen_plot.png"
plot2="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_Pressure_plot.png"

f= crystal_energies_per_volume
df=np.diff(f)
n= densities
dn=np.diff(n)
pressure = ((df/dn)*(n[0:num_rd]+dn/2))-(f[0:num_rd]+df/2)

# Plot Crystal Free Energy per sphere vs Reduced Density
plt.plot(data[:,0], data[:,1])
plt.title('Crystal Free Energy per sphere vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')
#plt.savefig('plotcopy1.png')  DEL/REF
plt.savefig(plot1)

plt.figure()

# Plot Pressure vs Reduced Density
plt.plot(n[0:num_rd]+dn/2, pressure)
plt.title('Pressure vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Pressure')
plt.savefig(plot2)

plt.show()

