#!/usr/bin/python3

# This program runs figs/new-melting.cpp for many different densities for 
# input parameters specified within this program.
# It is used to reproduce data, which is why the input paratmeters are 
# specified within the program, and are not intended to be changed.
# (Use compute-isotherm.py to run figs/new-melting.cpp for many different 
# densities and temperatures that are not fixed to reproduce published data.)
# ie. figs/new-melting.mkdat --kT 3 --n 1.39 --gwstart 0.01 --gwend 0.2 --gwstep 0.01 --fv 0 --dx 0.5 --seed 1 --mc-error 0.001 --mc-constant 5 --mc-prefactor 50000 --d data/phase-diagram --filename nm-kT_3-n_1.39.dat

# NOTE: Run this script from deft/papers/fuzzy-fmt
#       with command ./figs/compute-phasediagram.py

import numpy as np
import os
import argparse
import array as arr

def run_new_melting(kT, n, gwstart, gwend, gwstep, fv=0, dx=0.5, seed=5,
                    mcerror=1e-3, mcconstant=5, mcprefactor=50000,
                    avoid_rq=False):
    name = 'nm-kT_%g-n_%g_seed5' % (kT, n)
    cmd = 'rq run -J %s' % name
    if avoid_rq:
        cmd = ''
    cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gwstart %g --gwend %g --gwstep %g' % (gwstart, gwend, gwstep)
    cmd += ' --fv %g --dx %g --seed %g' % (fv, dx, seed)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    cmd += ' --d data/phase-diagram'
    cmd += ' --filename %s.dat' % name
    print(cmd)
    assert(os.system(cmd) == 0)

#kTs = np.arange(3, 0.05, -0.1)  #for small plot in paper
#kTs = np.arange(2.7, 0.05, -0.1)  #finish....then delete
kTs = np.arange(3, 0.05, -2.5)  #finish....then delete
#kTs = np.arange(0.5, 0.45, -0.1)  #trouble data points
densities = np.arange(0.01, 1.53, 0.02)

for kT in kTs:
    for n in densities:
        # It seems that for low temperature and high density we need to run a more
        # accurate and precise Monte Carlo, otherwise we end up seeing n_3 > 1 and get
        # NaNs.
        #if 0.499 < kT < 0.501 and 0.94 < (1/n) < 1.15:
        #if 0.499 < kT < 0.501 and 0.89 < (1/n) < 1.3:   #for small plot in paper
        if 0.499 < kT < 0.501 and 0.89 < (1/n) < 1:    #trouble data points
	    # We want high accuracy for this one!
            run_new_melting(kT, n, gwstart=0.01, gwend=0.11, gwstep=0.01, mcerror=1e-4) #trouble data points
        else:
            #pass # for now, skip these simulations that we already ran.
            run_new_melting(kT, n, gwstart=0.01, gwend=0.2, gwstep=0.01)

