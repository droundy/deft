#!/usr/bin/python3

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

parser = argparse.ArgumentParser(description='Creates data for plots of gw vs n and FE vs n.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--nmin', type=float,
                    help='min density - REQUIRED')
parser.add_argument('--nmax', type=float,
                    help='max density - REQUIRED')
parser.add_argument('--dn', type=float,
                    help='change of density', default=0.01)

parser.add_argument('--dgw', type=float,
                    help='change in gw', default=0.01)
parser.add_argument('--maxgw', type=float,
                    help='max gw', default=0.2)
parser.add_argument('--mingw', type=float,
                    help='min gw', default=0.01)

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
                    
parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')

args=parser.parse_args()

kT=args.kT

if args.dn:
    dn=args.dn
else :
    dn=0.01

if args.dgw:
    dgw=args.dgw
else :
    dgw=0.01

if args.maxgw:
    maxgw=args.maxgw
else :
    maxgw=0.2  
    
if args.mingw:
    mingw=args.mingw
else :
    mingw=0.01

if args.fv:
    fv=args.fv
else :
    fv=0

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
    
if args.tensor:    
    for n in np.arange(args.nmin, args.nmax, dn):
        cmd = 'rq run -J isotherm-kT%g-n%g-tensor' % (kT, n)
        cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
        cmd += ' --gwstart %g --gwend %g --gwstep %g' % (mingw, maxgw, dgw)
        cmd += ' --fv %g --dx %g' % (fv, dx)
        cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
        cmd += ' --filename isotherm-kT-%g_tensor.dat' % kT
        cmd += ' --tensor'
        print(cmd)
        os.system(cmd)
else :
    for n in np.arange(args.nmin, args.nmax, dn):
        cmd = 'rq run -J isotherm-kT%g-n%g' % (kT, n)
        cmd += ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
        cmd += ' --gwstart %g --gwend %g --gwstep %g' % (mingw, maxgw, dgw)
        cmd += ' --fv %g --dx %g' % (fv, dx)
        cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
        cmd += ' --filename isotherm-kT-%g.dat' % kT
        print(cmd)
        os.system(cmd)
