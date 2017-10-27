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
parser.add_argument('--d', metavar='     directory', type=str, default="crystalization",
                    help='directory for data files')                    

args=parser.parse_args()
print
print args  #for debug

densities=args.n
temperatures=args.t
nstart=args.nstart
nend=args.nend
nstep=args.nstep
tstart=args.tstart
tend=args.tend
tstep=args.tstep
data_dir=args.d
    
if args.nstart:
    densities = np.arange(nstart, nend+nstep, nstep, float)
    
if args.tstart:
    temperatures = np.arange(tstart, tend+tstep, tstep, float)

#If make a higher level script, move this question to the top level!
if args.d:
    print
    print "Data directory is deft/papers/fuzzy-fmt/"+data_dir
    print 
else:
    print
    print "Do you want to save default directory [fuzzy-fmt]/crystalization before it is over-written?"
    wait = raw_input("If not, press the ENTER key to continue program...")
    print 
    
print "Running new-melting with temperatures:", temperatures
print "and densities:", densities
print

for i in range(0,len(temperatures)):
    for j in range(0,len(densities)):
        print
        print "Temperature:", temperatures[i], "Density:", densities[j]  #testing for loop 
        os.system('figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10 --dir %s' %(temperatures[i],densities[j], data_dir))
        #os.system('rq run -J new-melting-%g-%g figs/new-melting.mkdat --kT %g --rd %g --fvstart 0.0 --fvend 1.0 --fvstep 0.2 --gwstart 0.01 --gwstep 10'
                  %(temperatures[i], densities[j], temperatures[i],densities[j]))

##----------------------------------------------------------------------
##NOTE: lattice_constant will be divided by gwstep     
   
        




