#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#if len(sys.argv) < 5:
#    print "Usage: T dfv dgw [DENSITIES TO COMPUTE]..."
#    exit(1)

#temp=sys.argv[1]    
#fvstep=sys.argv[2]
#gwstep=sys.argv[3]   #lattice_constant will be divided by this number
#rdensities=sys.argv[4:]  #if enter by command line and not by hand
#print rdensities
#temp_num=float(temp)

#os.system('rm newmeltdataout.dat')  //delete
#os.system('rm plot.dat')  //delete

#rdensities=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]  #if enter by hand and not by command line
rdensities=[0.2, 0.1]  #if enter by hand and not by command line TEST
num_rd=len(rdensities)-1 


#everything_file="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_dat.dat"  // delete
#best_file="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"  // delete from this file ... maybe include in plot script?

densities = np.array([float(d) for d in rdensities])
crystal_energiesdensities = np.zeros_like(densities)  #initializing...
crystal_energies_per_volume = np.zeros_like(densities)
energy_differences = np.zeros_like(densities)

for i in range(len(densities)):
    #if enter by hand and not by command line...
    rdensity=densities[i]  
    #print 'running with', rdensity     
    os.system('figs/new-melting.mkdat --kT 2 --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (rdensity))    
    #if enter by command line and not by hand...
#    rdensity_num=densities[i]   
#    print 'running with', rdensity_num 
#    os.system('figs/new-melting.mkdat --kT 2 --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % rdensity_num)  


#    os.system('cat *-bestdat.dat >> plot.dat')  # FIX THIS!

#    thisdata = np.loadtxt('kT%grd%g-bestdat.dat' % temp  % rdensity)
#    crystal_energies[i] = thisdata[1]
#    energy_differences[i] = thisdata[2]
#    crystal_energies_per_volume[i] = thisdata[3]

#print 'densities', densities
#print 'crystal_energies', crystal_energies
#print 'energy_differences', energy_differences
#print 'crystal_energies_per_volume', crystal_energies_per_volume

#cmdcpdat="cp %s/%s-alldat.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_dat.dat" 
#os.system(cmdcpdat)
#cmdcpplotdat="cp %s/%s-plot.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"  
#os.system(cmdcpplotdat)


#Plotting details...

#plot1="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_FEvsDen_plot.png"
#plot2="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_Pressure_plot.png"

#f= crystal_energies_per_volume
#df=np.diff(f)
#n= densities
#dn=np.diff(n)
#pressure = ((df/dn)*(n[0:num_rd]+dn/2))-(f[0:num_rd]+df/2)

# Plot Crystal Free Energy per sphere vs Reduced Density
#plt.plot(data[:,0], data[:,1])
#plt.title('Crystal Free Energy per sphere vs Reduced Density')
#plt.xlabel('Reduced Density')
#plt.ylabel('Free Energy')
#plt.savefig('plotcopy1.png')  DEL/REF
#plt.savefig(plot1)

#plt.figure()

# Plot Pressure vs Reduced Density
#plt.plot(n[0:num_rd]+dn/2, pressure)
#plt.title('Pressure vs Reduced Density')
#plt.xlabel('Reduced Density')
#plt.ylabel('Pressure')
#plt.savefig(plot2)

#plt.show()

