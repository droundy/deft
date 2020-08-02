#!/usr/bin/python2

#This program runs figs/new-melting.cpp with a particular temperature and density for many seeds

#NOTE: Run this script from directory deft/papers/fuzzy-fmt with the command ./run-10seeds.py 

#DOESN'T WORK THE WAY I WANTED IT TO - PROBLEM WITH FILENAMES!

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Creates data for error analysis.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
                    
parser.add_argument('--n', type=float,
                    help='density - REQUIRED')   
                                    
parser.add_argument('--fv', metavar='vacancies', type=float,
                    help='fraction of vacancies - Default 0')                   

parser.add_argument('--seeds', type=int,
                    help='number of seeds to run - REQUIRED')

parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')
args=parser.parse_args()

if args.fv:
    fv=args.fv
else :
    fv=0

def run_new_melting(kT, n, gwstart, gwend, seed, gwstep=0.01, fv=0, dx=0.5,
                    mcerror=1e-3, mcconstant=5, mcprefactor=50000,
                    avoid_rq=True):
    name = 'nm-kT_%g-n_%g_seed_%g' % (kT, n, seed)
    if args.tensor:
        name = 'tensor-'+name
    cmd = 'rq run -J %s' % name
    if avoid_rq:
		cmd = ''
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --d data/manyseeds'
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g --seed %g' % (fv, dx, seed)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    if args.tensor:
        cmd += ' --tensor'
    cmd += ' --filename %s.dat' % name
    print(cmd)
    assert(os.system(cmd) == 0)

if args.tensor:
    print("need to do this later")
    exit(0)

    
for seed in range(1, args.seeds+1):   
  print "seed=%g" % (seed)   
  run_new_melting(0.5, 0.96, 0.01, 0.2, seed)
       






