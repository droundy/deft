#!/usr/bin/python2

import os
import argparse
import numpy as np
import sys

#parser = argparse.ArgumentParser(description='Variables for running mew-melting:', epilog="stuff...")
parser = argparse.ArgumentParser(description="Must enter: EITHER a list of one or more denisties with --n OR density loop variables --nstart, --nend, --nstep; --t=2 by default otherwise enter EITHER a list of one or more temperatures with --t OR --tstart, --tend, --tstep.")
group1 = parser.add_mutually_exclusive_group(required=True)
group2 = parser.add_mutually_exclusive_group()
group1.add_argument('--n', metavar='density', type=float, nargs='+', 
                    help='reduced density(s)')
group2.add_argument('--t', metavar='temperature', type=float, nargs='+', default=[2],
                    help='reduced temperature(s) (default: 2)')
group1.add_argument('--nstart', metavar='nloop_start', type=float,
                    help='starting reduced density')
parser.add_argument('--nend', metavar='  nloop_end', type=float,
                    help='ending reduced density')
parser.add_argument('--nstep', metavar=' nloop_step', type=float, default=0.1,
                    help='step reduced density by (default: 0.1)')
group2.add_argument('--tstart', metavar='tloop_start', type=float,
                    help='starting temperature kT')
parser.add_argument('--tend', metavar='  tloop_end', type=float,
                    help='ending temperature kT')
parser.add_argument('--tstep', metavar=' tloop_step', type=float, default=1,
                    help='step temperature kT by (default: 1.0)')

args=parser.parse_args()
print args

densities=args.n
temperatures=args.t
nstart=args.nstart
nend=args.nend
nstep=args.nstep
tstart=args.tstart
tend=args.tend
tstep=args.tstep
    
if args.nstart:
    densities = np.arange(nstart, nend+nstep, nstep, float)
    
if args.tstart:
    temperatures = np.arange(tstart, tend+tstep, tstep, float)

for i in range(0,len(temperatures)):
    print temperatures[i]
    for j in range(0,len(densities)):
        print densities[j]
        #os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 
        #  --fvstep 0.2 --gwstart 0.01 --gwstep 10' %(temperatures[i],denstities[j])) 
        #print temperatures[i], densities[j]  #testing for loop
        
print "Running with densities: ", densities
print "length of densities=", len(densities)
print "Running with temperatures: ", temperatures
print "length of temperatures=", len(temperatures)
#print "nstart is", nstart 
#print "nend is", nend 
#print "nstep is", nstep
#print "tstart is", tstart
#print "tend is", tend 
#print "tstep is", tstep     

print "DONE"

##----------------------------------------------------------------------

##print "Do you want to save default directory deft/papers/fuzzy-fmt/crystalization 
#        before it is over-written?"
##os.system('rm ')  #ASK-remove data directory? - might not want to do this!

##NOTE: lattice_constant will be divided by gwstep     
   
        




