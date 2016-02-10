#!/usr/bin/python2
import numpy as np
import math
import string
import os

def dos_energy_ns(dbase, i):
    # read in files and compute U, F, and S given a directory and recursion level
    ln_dos = {}
    energy = {}
    ln_dos_hist = [] # need to initialize this correctly; how many elements?
    for N in range(2,100): #maybe require N as an argument?
        fbase = '%s/%01d/N%03d/lv-data-dos.dat' % (dbase, i, N) 
        if os.path.isfile(fbase+"-dos.dat"):
            ln_dos_hist = np.loadtext(fbase, ndmin=2)
            ln_dos[N] = ln_dos_hist[:,1]
            energy[N] = ln_dos_hist[:,0] #slicing may be wrong, need to look at this harder. Also, if a file doesn't exist, index is no longer =N...
        else:
            pass
    return ln_dos, energy

def F_hardsphere(dbase, i):
    fbase = dbase + '/absolute/'
    Fs = []  # initialize correctly or else use append....
    valid = [] # maybe use a list or class... these aren't square, so arrays seem a bad iead (lots of zeros) 
    failed = [] # maybe these don't even need to be stored?
    ratios = []
    Ns = []
    for j in range(2, 200):
        fname = fbase + '%05d' % (j)
        if os.path.isfile(fname+'.dat'):
            # if file exists, load the text from the file
            pass
        for k in range(0,len(valid)):
            # compute absolute Fs using read in data. Fs indexed by [N]
            ratios[k,j] = failed[k,j] / valid[k,j]   # indexed by [filenumber, N]
            Fs[j] += -np.log(ratios[k,j])/ Ns[j]
            
def U_F_S_ns(dbase, i):
    ln_dos, energy = dos_energy_ns(dbase, i)
    Ns = eval(ln_dos.keys())
    T_bins, dT, T_range = 1e3, 1/T_bins, np.arange(dT, 1, dT) # change Tmax # Generate array of temperatures
    Z = np.zeros(len(T_range), len(Ns))
    Zinf, U, F, S = (np.zeros_like(Z) for i in range(4))
    
    for N in set(Ns):
        j=0
        for T in set(T_range):
            k=0
            ln_dos_boltz = ln_dos[N] - energy[N]/T
            dos_boltz = np/exp(ln_dos_boltz - ln_dos_boltz.max())
            Z[k,j] = sum(dos_boltz)  #indexed by [T, N]
            
            U[k,j] = sum[energy[N]*dos_boltz]/Z[k,j]
            F[k,j] = # need to map the [:,N] elements with absolute calculations... also 
            
    
    
    
            
    
