#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#for arg in sys.argv[1:]:
#    print arg

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

for rdensity in rdensities:
    rdensity_num=float(rdensity)   #if enter by command line and not by hand
#   os.system('figs/new-melting.mkdat 2 %g -1 -1' % (rdensity))    #if enter by hand and not by command line
    os.system('figs/new-melting.mkdat 2 %g -1 -1' % rdensity_num)  #if enter by command line and not by hand
    os.system('cat newmeltbestdata.dat >> plot.dat')


cmdcpdat="cp newmeltdataout.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_newmeltdataout.dat" 
os.system(cmdcpdat)
cmdcpplotdat="cp plot.dat DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_plot.dat"  
os.system(cmdcpplotdat)

plot1="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_FEvsDen_plot.png"
plot2="DATA_KIR/NewStuff/gw"+gwstep+"fv"+fvstep+"_Pressure_plot.png"

data = np.loadtxt('plot.dat')
f=data[:,3]
df=np.diff(data[:,3])
n=data[:,0]
dn=np.diff(data[:,0]) 
pressure = ((df/dn)*(n[0:num_rd]+dn/2))-(f[0:num_rd]+df/2)

# Plot Crystal Free Energy per sphere vs Reduced Density
plt.plot(data[:,0], data[:,1])
plt.title('Crystal Free Energy per sphere vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Free Energy')
#plt.savefig('plotcopy1.png')  DEL/REF
plt.savefig(plot1)
plt.show()

# Plot Pressure vs Reduced Density
plt.plot(n[0:num_rd]+dn/2, pressure)
plt.title('Pressure vs Reduced Density')
plt.xlabel('Reduced Density')
plt.ylabel('Pressure')
plt.savefig(plot2)
plt.show()

