#!/usr/bin/python3

#This program runs figs/new-melting.cpp for many different densities for input parameters specified within this program
#It is used to reproduce data, which is why the input paratmeters are specified within the program, and are not intended to be changed.
#(use compute-isotherm.py to run figs/new-melting.cpp for many different densities based on user input)
#ie. figs/new-melting.mkdat --kT 0.5 --n 1.06 --gwstart 0.01 --gwend 0.2 --gwstep 0.01 --fv 0.01 --dx 0.5 --mc-error 0.001 --mc-constant 5 --mc-prefactor 50000 --filename isotherm-kT-0.5_tensor.dat --tensor

#NOTE: Run this script from directory deft/papers/fuzzy-fmt with command python3 ./compute-phasediagram.py --tensor(optional)
#Currently, data is stored in the direction deft/papers/fuzzy-fmt/newdata_tensor/phase-diagram5



import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse
import array as arr

parser = argparse.ArgumentParser(description='Creates data for phase diagrams.')

parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')
args=parser.parse_args()

def run_new_melting(kT, n, gwstart, gwend, gwstep, fv=0, dx=0.5, seed=1,
                    mcerror=1e-3, mcconstant=5, mcprefactor=50000,
                    avoid_rq=False):
    name = 'nm-kT_%g-n_%g' % (kT, n)
    if args.tensor:
        name = 'tensor-'+name
    cmd = 'rq run -J %s' % name
    if avoid_rq:
        cmd = ''
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g --seed %g' % (fv, dx, seed)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    if args.tensor:
        cmd += ' --d newdata_tensor/phase-diagram5'   #once new results are compared to the old, remove the "2" at the end of phase-diagram
        cmd += ' --tensor'
    else:
        #cmd += ' --d data/phase-diagram'
        cmd += ' --d newdata/phase-diagram5' #once new results are compared to the old, remove the "2" at the end of phase-diagram
    cmd += ' --filename %s.dat' % name
    print(cmd)
    assert(os.system(cmd) == 0)

#kTs = np.arange(2, 0.3, -0.1)   # kT less than 0.4 don't have solutions
#kTs = np.arange(20, 2, -2) 
kTs = np.arange(.5, .1, -.5) 
#kTs = np.arange(200, 20, -20)
#kTs = np.arange(.05, .005, -0.02)
#kTs = np.arange(3, 0.05, -0.1)
#kTs = np.arange(200, 100, -20)
#kTs = np.arange(1.15, 0.7, -0.05)
#kTs=np.append(kTs, 0.01)
for kT in kTs:
    #for n in np.arange(0.01, 0.59, 0.02):   #homogeneous fluid
        #run_new_melting(kT, n, 0.001, 0.001, 0.001, avoid_rq=True)    #temp, density, gw_start, gw_end, gw_step
    # for n in np.arange(0.59, 1.2, 0.02):    #crystal
        # if args.tensor:
            # run_new_melting(kT, n, 0.01, 0.2, 0.01)    #temp, density, gw_start, gw_end, gw_step
        # else:
            # run_new_melting(kT, n, 0.01, 0.2, 0.01)
    for n in np.arange(.6, 1.4, 0.2):
         #run_new_melting(kT, n, 0.001, 0.001, 0.001) #super fast for homogeneous at the high temps I'm running
         run_new_melting(kT, n, 0.01, 0.2, 0.01)
    # for n in np.arange(1.22, 1.42, 0.02):
        # run_new_melting(kT, n, 0.01, 0.2, 0.01)
    #for n in np.arange(1.42, 1.82, 0.02):
        #run_new_melting(kT, n, 0.01, 0.2, 0.01)
    #for n in np.arange(1.82, 2.8, 0.2):
        #run_new_melting(kT, n, 0.01, 0.2, 0.01)

