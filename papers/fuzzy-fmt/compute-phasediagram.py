#!/usr/bin/python3

#This program runs figs/new-melting.cpp for many different densities for input parameters specified within this program
#It is used to reproduce data, which is why the input paratmeters are specified within the program, and are not intended to be changed.
#(use compute-isotherm.py to run figs/new-melting.cpp for many different densities based on user input)
#ie. figs/new-melting.mkdat --kT 0.5 --n 1.06 --gwstart 0.01 --gwend 0.2 --gwstep 0.01 --fv 0.01 --dx 0.5 --mc-error 0.001 --mc-constant 5 --mc-prefactor 50000 --d newdata_tensor/phase-diagram4 --filename isotherm-kT-0.5_tensor.dat

#NOTE: Run this script from directory deft/papers/fuzzy-fmt with command ./compute-phasediagram.py
#Currently, data is stored in the directory deft/papers/fuzzy-fmt/newdata_tensor/phase-diagram4

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse
import array as arr

def run_new_melting(kT, n, gwstart, gwend, gwstep, fv=0, dx=0.5, seed=1,
                    mcerror=1e-3, mcconstant=5, mcprefactor=50000,
                    avoid_rq=False):
    name = 'nm-kT_%g-n_%g' % (kT, n)
    cmd = 'rq run -J %s' % name
    if avoid_rq:
        cmd = ''
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g --seed %g' % (fv, dx, seed)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    cmd += ' --d newdata_tensor/phase-diagram4' 
    cmd += ' --filename %s.dat' % name
    print(cmd)
    assert(os.system(cmd) == 0)

#kTs = np.arange(2, 0.3, -0.1)   # kT less than 0.4 don't have solutions
#kTs = np.arange(20, 2, -2) 
#kTs = np.arange(39, 0.5, -2) 
#kTs = np.arange(200, 20, -20)
#kTs = np.arange(.05, .005, -0.02)
###kTs = np.arange(3, 0.05, -0.1)  #for small plot in paper
kTs = np.arange(40, 0.05, -40)  #for small plot in paper
#kTs = np.arange(200, 100, -20)
#kTs = np.arange(1.15, 0.7, -0.05)
#kTs=np.append(kTs, 0.01)
for kT in kTs:
    ###for n in np.arange(0.01, 0.59, 0.02):   #homogeneous fluid                            #for small plot in paper
        #run_new_melting(kT, n, gwstart=0.001, gwend=0.001, gwstep=0.001, avoid_rq=True)  #super fast for homogeneous free energy (not used for paper)
        ###run_new_melting(kT, n, gwstart=0.01, gwend=0.2, gwstep=0.01)                      #for small plot in paper
    ###for n in np.arange(0.59, 1.2, 0.02):    #crystal                                      #for small plot in paper
        ###run_new_melting(kT, n, gwstart=0.01, gwend=0.2, gwstep=0.01)                      #for small plot in paper
    ##for n in np.arange(0.60, 2.32, 0.02):
         #run_new_melting(kT, n, 0.001, 0.001, 0.001) #super fast for homogeneous at the high temps I'm running
         ##run_new_melting(kT, n, gwstart=0.01, gwend=0.2, gwstep=0.01)
    # for n in np.arange(1.22, 1.42, 0.02):
        # run_new_melting(kT, n, 0.01, 0.2, 0.01)
    #for n in np.arange(1.42, 1.82, 0.02):
        #run_new_melting(kT, n, 0.01, 0.2, 0.01)
    #for n in np.arange(1.82, 2.8, 0.2):
        #run_new_melting(kT, n, 0.01, 0.2, 0.01)
    for n in np.arange(0.01, 1.2, 0.02):    #crystal                                      #for small plot in paper
        run_new_melting(kT, n, gwstart=0.01, gwend=0.2, gwstep=0.01)  

