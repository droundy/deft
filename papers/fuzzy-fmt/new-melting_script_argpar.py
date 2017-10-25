#!/usr/bin/python2

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys



parser = argparse.ArgumentParser(description='Variables for running mew-melting.')
parser.add_argument('--n', metavar='n', type=float, nargs='+',
                    help='reduced density(s) (default: 1.2)')
parser.add_argument('--t', metavar='t', type=float, 
                    help='temperature (default: 2)')
parser.add_argument('--nstart', metavar='n_start', type=float, 
                    help='starting reduced density (default: 0.2)')
parser.add_argument('--nend', metavar='n_end', type=float, 
                    help='ending reduced density (default: 1.8)')
parser.add_argument('--nstep', metavar='n_step', type=float, 
                    help='step reduced density by (default: 0.1)')
parser.add_argument('--tstart', metavar='t_start', type=float, 
                    help='starting temperature kT (default: 2.0)')
parser.add_argument('--tend', metavar='t_end', type=float, 
                    help='ending temperature (default: 2.0)')
parser.add_argument('--tstep', metavar='ntstep', type=float, 
                    help='step temperature by (default: 1.0)')

args=parser.parse_args()
#print args

n=args.n
print "n is", n 
temp=args.t
print "temp is", temp
nstart=args.nstart
print "nstart is", nstart
nend=args.nend
print "nend is", nend
nstep=args.nstep
print "nstep is", nstep
tstart=args.tstart
print "tstart is", tstart
tend=args.tend
print "tend is", tend
tstep=args.tstep
print "tstep is", tstep


print "DONE"


##-----------------------------
#if len(sys.argv) < 3:
#    print "Usage: T dfv dgw [DENSITIES TO COMPUTE]..."  #no longer needed?
#    exit(1)

##make these loops over T and rdensities?
#temp=sys.argv[1]    
#rdensities=sys.argv[2:]  #if enter by command line and not by hand
#print rdensities
#temp_num=float(temp)

##print "Do you want to save default directory deft/papers/fuzzy-fmt/crystalization before it is over-written?"
##os.system('rm ')  #ASK-remove data directory? - might not want to do this!

##rdensities=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]  #if enter by hand and not by command line
##rdensities=[0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]  #if enter by hand and not by command line
#num_rd=len(rdensities)-1 

##ASK-change rdensities strings from command line input to float numbers
#densities = np.array([float(d) for d in rdensities])
  
#for i in range(len(densities)):
   ##if enter by hand and not by command line...
    ##rdensity=densities[i]  
    ##print 'running with', rdensity     
    ##os.system('figs/new-melting.mkdat --kT 2 --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (rdensity))    
   ##if enter by command line and not by hand...
    #rdensity_num=densities[i]   
    #print 'running with', rdensity_num  
    #os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (temp_num, rdensity_num))
      

##NOTE: lattice_constant will be divided by gwstep     
   
        




