#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./nm_FE_vs_gw_generate_data.py [filename.dat] --kT [required] --n [required] --fv [] --dx [] --mcerror [] --mcconstant [] --mcprefactor [] --seeds []
#For help type: ./nm_FE_vs_gw_generate_data.py --help
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")

#parser.add_argument('filedat', metavar='data filename', type=str,
#                    help='name of datafile to store FE vs gw data') 
parser.add_argument('--kT', metavar='reduced temperature', type=float,
                    help='reduced temperature - required') 
parser.add_argument('--n', metavar='reduced density', type=float,
                    help='reduced density - required') 
parser.add_argument('-fv', metavar='fraction of vacancies', type=float,
                    help='fraction of vacancies')  
parser.add_argument('--dx', metavar='dx', type=float,
                    help='scaling dx') 
#parser.add_argument('--mcerror', metavar='mc_error', type=float,
#                    help='mc_error') 
parser.add_argument('--mcconstant', metavar='mc_constant', type=int,
                    help='mc_constant') 
parser.add_argument('--mcprefactor', metavar='mc_prefactor', type=int,
                    help='mc_prefactor')
parser.add_argument('--seeds', metavar='seeds', type=int,
                    help='number of seeds ran') 
                                        
args=parser.parse_args()

kT=args.kT
n=args.n
#seeds=args.seeds

if args.fv:
    fv=args.fv
else :
    fv=0

if args.seeds:
    seeds=args.seeds
else :
    seeds=50

if args.dx:
    dx=args.dx
else :
    dx=.5

#if args.mcerror:
#    mcerror=args.mcerror
#else :
#    mcerror=0.001
#fix formatting in data file out name!

if args.mcconstant:
    mcconstant=args.mcconstant
else :
    mcconstant=5
    
if args.mcprefactor:
    mcprefactor=args.mcprefactor
else :
    mcprefactor=50000
    
#if args.filedat:    
#    datafile=args.filedat+'.dat'

for gw in [.4, .5, .6, .7, .8, .9]:
#for gw in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3]: 
  print "gw=%g" % (gw)
  for seed in range(1, seeds+1):   
    print "seed=%g" % (seed)
    os.system('figs/new-melting.mkdat --kT 2 --n 1.3 --gw %g  --fv 0 --dx %g  --mc-error .001   --mc-constant %g --mc-prefactor %g --seed %g | tail -n 2   >> FE_vs_gw_kT%g_n%g_fv0_dx%g_mcerror0.001_mcconstant%g_mcprefactor%g.dat' % (gw, dx, mcconstant, mcprefactor, seed, kT, n, dx, mcconstant, mcprefactor)) 

    

