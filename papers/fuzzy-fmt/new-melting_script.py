#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./new-melting_script.py followed by agruments as 
#described with --usage flag.

import os
import argparse
import numpy as np
import sys

#parser = argparse.ArgumentParser(description='Variables for running mew-melting:', epilog="stuff...")
parser = argparse.ArgumentParser(description="Must enter: EITHER a list of one or more denisties with --n OR density loop variables --nstart, --nend, --nstep; --t=2 by default otherwise enter EITHER a list of one or more temperatures with --t OR --tstart, --tend, --tstep.")
group1 = parser.add_mutually_exclusive_group(required=True)
group2 = parser.add_mutually_exclusive_group()
group3 = parser.add_mutually_exclusive_group()
group4 = parser.add_mutually_exclusive_group()
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
parser.add_argument('--d', metavar='     directory', type=str, default="crystallization",
                    help='directory for data files') 
parser.add_argument('--fv', metavar=' loop option', type=float, default=-1,   
                    help='gwidth or -1 loop (default: -1)')
parser.add_argument('--fvstart', metavar='  fvloop_start', type=float, default=0,   
                    help='starting fv (default: 0)')
parser.add_argument('--fvend', metavar='  fvloop_end', type=float, default=1,
                    help='ending gwidth (default: 1)')
parser.add_argument('--fvstep', metavar='fvloop_step', type=float, default=0.2,
                    help='step fv by (default: 0.2)') 
parser.add_argument('--gw', metavar=' loop option', type=float, default=-1,   
                    help='gwidth or -1 loop without lattice ref or -2 loop with lattice ref (default: -1)')
parser.add_argument('--gwstart', metavar='  gwloop_start', type=float, default=0.01,   
                    help='starting gwidth (default: 0.01)')
group3.add_argument('--gwend', metavar='  gwloop_end', type=float,
                    help='ending gwidth')
group4.add_argument('--gwstep', metavar=' gwloop_step', type=float,
                    help='step gwidth') 
group3.add_argument('--gwlend', metavar=' gwloop_latend', type=float, default=0.5,
                    help='ending gwidth will be computed lattice_constant multiplied by this number (default: 0.5)')
group4.add_argument('--gwlstep', metavar='gwloop_latstep', type=float, default=0.1,
                    help='gwidth will step by computed lattice_constant multiplied by this number (default: 0.1)') 
parser.add_argument('--dx', metavar='dx grid spacing', type=float, default=0.01,
                    help='grid spacing (default: 0.01)')                   
                                   

args=parser.parse_args()
print
print args  #also put this in data record of what ran!

densities=args.n
temperatures=args.t
nstart=args.nstart
nend=args.nend
nstep=args.nstep
tstart=args.tstart
tend=args.tend
tstep=args.tstep
data_dir=args.d
dx=args.dx

fv_start=args.fvstart
fv_end=args.fvend
fv_step=args.fvstep
gw=args.gw
gwidth_start=args.gwstart
gwidth_end=args.gwend
gwidth_step=args.gwstep
gwidth_latend=args.gwlend
gwidth_latstep=args.gwlstep
    
if args.nstart:
    densities = np.arange(nstart, nend, nstep, float)
    
if args.tstart:
    temperatures = np.arange(tstart, tend, tstep, float)
    
print     
print "Running new-melting_script.py with:"
print "  kT Temperatures:", temperatures
print "  Reduced Densities:", densities
print "  Data directory: deft/papers/fuzzy-fmt/"+data_dir
print 


for i in range(0,len(temperatures)):
    for j in range(0,len(densities)):
        print
        print "Temperature:", temperatures[i], "Density:", densities[j]
        cmd = ''
        if os.system('which rq') == 0:
            cmd = 'rq run -J new-melting-kT=%g-n=%g ' % (temperatures[i],densities[j])
        cmd += 'figs/new-melting.mkdat --kT %g --n %g' % (temperatures[i],densities[j])
        cmd += ' --d %s --dx' % (data_dir, dx)
        cmd += ' --fvstart %g --fvend %g --fvstep %g --fv %g' % (fv_start, fv_end, fv_step, fv)
        if args.gwend or args.gwstep:
            cmd += ' --gwstart %g --gwend %g --gwstep %g --gw %g' % (gwidth_start, gwidth_end, gwidth_step, gw)
        else: 
            cmd += ' --gwstart %g --gwlend %g --gwlstep %g' % (gwidth_start, gwidth_latend, gwidth_latstep)

        print(cmd)
        os.system(cmd)

   
   
        




