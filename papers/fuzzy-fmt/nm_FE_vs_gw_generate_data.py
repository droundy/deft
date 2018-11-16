#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./nm_FE_vs_gw_generate_data.py --kT [required] --n [required] --fv [] --dx [] --mcerror [] --mcconstant [] --mcprefactor [] --seeds []  --datafile []
#For help type: ./nm_FE_vs_gw_generate_data.py --help


import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Creates data for a FE vs gw plot.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED') 
parser.add_argument('--n', metavar='density', type=float,
                    help='reduced density - REQUIRED') 
                    
parser.add_argument('--fv', metavar='vacancies', type=float,
                    help='fraction of vacancies - Default 0')  
parser.add_argument('--gw', metavar='width', type=float,
                    help='width of Gaussian - Default 0.01')  
parser.add_argument('--dx', metavar='dx', type=float,
                    help='scaling dx - Default 0.5') 
parser.add_argument('--mcerror', metavar='mc_error', type=float,
                    help='monte carlo mc_error - Default 0.001') 
parser.add_argument('--mcconstant', metavar='const', type=int,
                    help='monte carlo integration mc_constant - Default 5') 
parser.add_argument('--mcprefactor', metavar='prefac', type=int,
                    help='monte carlo integration mc_prefactor - Default 50000')
parser.add_argument('--seeds', metavar='seeds', type=int,
                    help='number of seeds ran - Default 1') 
parser.add_argument('--datafile', metavar='addtoname', type=str,
                    help='added to name of data file - Default "FE_vs_gw"') 
                    
                                        
args=parser.parse_args()

kT=args.kT
n=args.n

if args.fv:
    fv=args.fv
else :
    fv=0
    
if args.gw:
    gw=args.gw
else :
    gw=0.01

if args.dx:
    dx=args.dx
else :
    dx=.5

if args.mcerror:
    mcerror=args.mcerror
else :
    mcerror=0.001

if args.mcconstant:
    mcconstant=args.mcconstant
else :
    mcconstant=5
    
if args.mcprefactor:
    mcprefactor=args.mcprefactor
else :
    mcprefactor=50000
    
if args.seeds:
    seeds=args.seeds
else :
    seeds=1
    
if args.datafile:    
    datafile='FE_vs_gw_'+args.filedat
else :
    datafile='FE_vs_gw'
 
#for gw in [0.05]: 
#for gw in [.0000001, .000001, .00001, .0001, .001, .004, .006, .008]:
for gw in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]: #, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5]: 
#for gw in [.6, .7, .8]:
  print "gw=%g" % (gw)
  for seed in range(1, seeds+1):   
    print "seed=%g" % (seed)
    # os.system('rq run -J gw-%g-seed-%d figs/new-melting.mkdat --kT %g --n %g --gw %g  --fv %g --dx %g  --mc-error %g   --mc-constant %g --mc-prefactor %g --seed %g --filename %s_kT%g_n%g_fv%g_dx%g_mcerror%g_mcconstant%g_mcprefactor%g_seeds%g.dat'
    #           % (gw, seeds,
    #              kT, n, gw, fv, dx, mcerror, mcconstant, mcprefactor, seed, datafile, kT, n,
    #              fv, dx, mcerror, mcconstant, mcprefactor, seeds))
    #roundy#cmd = ('figs/new-melting.mkdat --kT %g --n %g --gw %g  --fv %g --dx %g  --mc-error %g   --mc-constant %g --mc-prefactor %g --seed %g --filename %s_kT%g_n%g_fv%g_dx%g_mcerror%g_mcconstant%g_mcprefactor%g_seeds%g.dat &'
    cmd = ('figs/new-melting.mkdat --kT %g --n %g --gw %g  --fv %g --dx %g  --mc-error %g   --mc-constant %g --mc-prefactor %g --seed %g --filename %s_kT%g_n%g_fv%g_dx%g_mcerror%g_mcconstant%g_mcprefactor%g_seeds%g.dat'          
              % (kT, n, gw, fv, dx, mcerror, mcconstant, mcprefactor, seed, datafile, kT, n,
                 fv, dx, mcerror, mcconstant, mcprefactor, seeds))

    print(cmd)
    os.system(cmd)



