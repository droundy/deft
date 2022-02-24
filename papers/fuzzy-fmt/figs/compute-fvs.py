#!/usr/bin/python3

# This program runs figs/new-melting.cpp for many different densities for 
# input parameters specified within this program.
# It is used to reproduce data, which is why the input paratmeters are 
# specified within the program, and are not intended to be changed.
# (Use compute-isotherm.py to run figs/new-melting.cpp for many different 
# densities and temperatures that are not fixed to reproduce published data.)
# ie. figs/new-melting.mkdat --kT 3 --n 1.39 --gwstart 0.01 --gwend 0.2 --gwstep 0.01 --fv 0 --dx 0.5 --seed 1 --mc-error 0.001 --mc-constant 5 --mc-prefactor 50000 --d data/phase-diagram --filename nm-kT_3-n_1.39.dat

# NOTE: Run this script from deft/papers/fuzzy-fmt
#       with command ./figs/compute-fvs.py  --kT [temp] --n [density]  --mcerror [1e-3 or 1e-4]
#ie.    ./figs/compute-fvs.py  --kT 0.5 --n 1.07  --mcerror 1e-3

import numpy as np
import os
import argparse
import array as arr

parser = argparse.ArgumentParser(description='Creates a plot of FEdiff vs gw specified temperature, density, and mcerror for a determined set of fraction of vacancies and seeds 1-5.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', metavar='density', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--mcerror', metavar='mcerror', type=float,
                    help='enter 1e-3 or 1e-4 - REQUIRED')

args=parser.parse_args()

kT=args.kT
n=args.n
mcerror=args.mcerror   

if mcerror == 1e-3:
   directory_name='fv-mcerr3'
if mcerror == 1e-4:
   directory_name='fv-mcerr4'               

def run_new_melting(kT, n, fv, seed, mcerror, gwstart, gwend, gwstep, dx=0.5, 
                    mcconstant=5, mcprefactor=50000,
                    avoid_rq=False):
    name = 'nm-kT_%g-n_%g_fv_%g_%s_seed%g' % (kT, n, fv, directory_name, seed)    
    cmd = 'rq run -J %s' % name
    if avoid_rq:
        cmd = ''
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g --seed %g' % (fv, dx, seed)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    cmd += ' --d data/phase-diagram-test-%s'  % (directory_name)
    cmd += ' --filename %s.dat' % name
    print(cmd)
    assert(os.system(cmd) == 0)
    
fvs = (0, 1e-9, 1e-8, 1e-6, 1e-4, 1e-2)
seeds = (1, 2, 3, 4, 5)

for seed in seeds:
 for fv in fvs:
     run_new_melting(kT, n, fv, seed, mcerror, gwstart=0.01, gwend=0.05, gwstep=0.04) 
     run_new_melting(kT, n, fv, seed, mcerror, gwstart=0.04, gwend=0.1, gwstep=0.04)  
           
