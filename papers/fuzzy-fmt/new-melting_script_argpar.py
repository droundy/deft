#!/usr/bin/python2

import os
import argparse
import numpy as np
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
print "density is", n
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

for i in range(0,len(n)):
    print "n[i] is", n[i]
   


##-----------------------------
#temp_num=float(temp)

##print "Do you want to save default directory deft/papers/fuzzy-fmt/crystalization before it is over-written?"
##os.system('rm ')  #ASK-remove data directory? - might not want to do this!

#num_rd=len(rdensities)-1 

##ASK-change rdensities strings from command line input to float numbers
#densities = np.array([float(d) for d in rdensities])
#if n=none
#    print (nend-nstart)/nstep
#    densities=
#    else
    #densities = np.array([float(d) for d in n])  #??

  
#for i in range(len(densities)):
    #rdensity_num=densities[i]   
    #print 'running with', rdensity_num 
#    print densities[i] 
    #os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (temp, rdensity_num)
    #os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (temp, densities[i]))   ??   

##NOTE: lattice_constant will be divided by gwstep     
   
        




