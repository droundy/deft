#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#if len(sys.argv) < 5:
#    print "Usage: T dfv dgw [DENSITIES TO COMPUTE]..."
#    exit(1)

#keep this for now...not sure what arguments we might want yet!
 #temp=sys.argv[1]    
 #fvstep=sys.argv[2]
 #gwstep=sys.argv[3]   #lattice_constant will be divided by this number
 #rdensities=sys.argv[4:]  #if enter by command line and not by hand
 #print rdensities
 #temp_num=float(temp)

#os.system('rm ')  #ASK-remove data directory? - might not want to do this!

#rdensities=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]  #if enter by hand and not by command line
rdensities=[0.3, 0.4]  #if enter by hand and not by command line TEST
num_rd=len(rdensities)-1 

#ASK-change rdensities strings from command line input to float numbers
densities = np.array([float(d) for d in rdensities])
  
for i in range(len(densities)):
   #if enter by hand and not by command line...
    rdensity=densities[i]  
    #print 'running with', rdensity     
    os.system('figs/new-melting.mkdat --kT 2 --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (rdensity))    
   #if enter by command line and not by hand...
    #rdensity_num=densities[i]   
    #print 'running with', rdensity_num 
    #os.system('figs/new-melting.mkdat --kT 2 --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % rdensity_num)  


#Do we need these in the plot file? - ASK!
 #crystal_energiesdensities = np.zeros_like(densities)  #initializing...
 #crystal_energies_per_volume = np.zeros_like(densities)
 #energy_differences = np.zeros_like(densities)
