#!/usr/bin/python2

#This program runs figs/new-melting.cpp for many different densities
#ie. figs/new-melting.mkdat --kT 0.5 --n 1.06 --gwstart 0.01 --gwend 0.2 --gwstep 0.01 --fv 0.01 --dx 0.5 --mc-error 0.001 --mc-constant 5 --mc-prefactor 50000 --filename isotherm-kT-0.5_tensor.dat --tensor


#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./compute-isotherm.py --kT [temp] --nmin [starting density] --nmax [ending density] --tensor(optional)
#For list of the many other options enter ./compute-isotherm.py --help


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Creates data for phase diagram.')

parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')
args=parser.parse_args()

def run_new_melting(kT, n, gwstart, gwend, gwstep=0.01, fv=0, dx=0.5, mcerror=1e-3, mcconstant=5, mcprefactor=50000):
    name = 'nm-kT_%g-n_%g' % (kT, n)
    if args.tensor:
        name = 'tensor-'+name
    cmd = 'rq run -J %s' % name
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' -d data/phase-diagram'
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g' % (fv, dx)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    if args.tensor:
        cmd += ' --tensor'
    cmd += ' --filename %s.dat' % name
    print(cmd)
    # os.system(cmd)

if args.tensor:
    print("need to do this later")
    exit(0)

for kT in np.arange(0.1, 1.15, 0.05):
    for n in np.arange(0.7, 1.0, 0.01):
        run_new_melting(kT, n, 0.01, 0.2)

    for n in np.arange(0.01, 0.69, 0.01):
        run_new_melting(kT, n, 0.001, 0.001)
