#!/usr/bin/python2

import os
import argparse
import numpy as np
import sys

parser = argparse.ArgumentParser(description='Variables for running mew-melting.')
parser.add_argument('--n', metavar='n', type=float, nargs='+', default=[1.2],
                    help='reduced density(s) (default: 1.2)')
parser.add_argument('--t', metavar='t', type=float, nargs='+', default=[2],
                    help='temperature(s) (default: 2)')
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
print args

densities=args.n
temperatures=args.t
nstart=args.nstart
nend=args.nend
nstep=args.nstep
tstart=args.tstart
tend=args.tend
tstep=args.tstep

if args.nstart and args.n:
    print "Error - EITHER enter --n OR --nstart, --nend, and --nstep"
if args.tstart and args.t:
    print "Error - EITHER enter --t OR --tstart, --tend, and --tstep"
    
if args.nstart:
    densities = np.arange(nstart, nend+nstep, nstep, float)
    
if args.tstart:
    temperatures = np.arange(tstart, tend+tstep, tstep, float)

for i in range(0,len(temperatures)):
    print temperatures[i]
    for j in range(0,len(densities)):
        print densities[j]
        #os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10' % (temperatures[i], denstities[j])) 
        
        
print "Running with densities: ", densities
#print "length of densities=", len(densities)
print "Running with temperatures: ", temperatures
#print "length of temperatures=", len(temperatures)
#print "nstart is", nstart 
#print "nend is", nend 
#print "nstep is", nstep
#print "tstart is", tstart
#print "tend is", tend 
#print "tstep is", tstep     

print "DONE"

##----------------------------------------------------------------------

##print "Do you want to save default directory deft/papers/fuzzy-fmt/crystalization before it is over-written?"
##os.system('rm ')  #ASK-remove data directory? - might not want to do this!

##NOTE: lattice_constant will be divided by gwstep     
   
        




