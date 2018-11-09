#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./nm_FE_vs_gw_generate_data.py [filename.dat] --kT [] --n [] --seeds []

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

#parser.add_argument('filedat', metavar='data filename', type=str,
#                    help='name of datafile to store FE vs gw data') 
parser.add_argument('--kT', metavar='reduced temperature', type=float,
                    help='reduced temperature') 
parser.add_argument('--n', metavar='reduced density', type=float,
                    help='reduced density') 
#parser.add_argument('--seeds', metavar='seeds', type=float,
#                    help='number of seeds ran') 
                                        
args=parser.parse_args()


kT=args.kT
n=args.n
#seeds=args.seeds
#datafile=args.filedat+'.dat'

seeds=1
mcprefactor=50000
mcconstant=5 
dx=.5

for gw in [.01, .02]: #test
#for gw in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]: 
  print "gw=%g" % (gw)
  for seed in range(1, seeds+1):   
    print "seed=%g" % (seed)
    os.system('figs/new-melting.mkdat --kT 2 --n 1.3 --gw %g  --fv 0 --dx %g  --mc-error .001   --mc-constant %g --mc-prefactor %g --seed %g | tail -n 2   >> FE_vs_gw_kT%g_n%g_dx%g_mcerror0.001.dat' % (gw, dx, mcconstant, mcprefactor, seed, kT, n, dx)) 
    #os.system('figs/new-melting.mkdat --kT %g --n %g --gw %g  --fv 0 --dx .5  --mc-error .001   --mc-constant %g --mc-prefactor %g --seed %g | tail -n 2   >> Hist_%g_%g_gw%g_mcerror0.001.dat' % (kT, n, gw, mcconstant, mcprefactor, seed, mcconstant, mcprefactor, gw)) 
    
    
    
  #os.system('./Histogram_original.py Hist_%g_%g_gw%g_mcerror0.001.dat --gw %g --mcconstant %g --seeds %g >> Hist_%g_mcerror0.001.dat' % (mcconstant, mcprefactor, gw, gw, mcconstant, seeds, mcprefactor)) 

