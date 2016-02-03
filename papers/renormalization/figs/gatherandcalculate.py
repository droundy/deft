#!/usr/bin/python2
import numpy as np
import math
import string
import os

def dos_energy_ns(dbase, i):
    # read in files and compute U, F, and S given a directory and recursion level
    ln_dos = {}
    energy = {}
    ln_dos_hist = []
    for N in range(2,100): #maybe require N as an argument?
        fbase = '%s/%01d/N%03d/lv-data-dos.dat' % (dbase, i, N) 
        if os.path.isfile(fbase+"-dos.dat"):
            ln_dos_hist = np.loadtext(fbase, ndmin=2)
            ln_dos[N] = ln_dos_hist[:,1]
            energy[N] = ln_dos_hist[:,0] #slicing may be wrong, need to look at this harder. Also, if a file doesn't exist, index is no longer =N...
        else:
            pass
    return ln_dos, energy

def dos_energy(fname):
    ln_dos = []
    energy = []
    
    ln_dos_hist = np.loadtext(fname, ndmin=2)
    ln_dos = ln_dos_hist[:,1]
    energy = ln_dos_hist[:,0]
    
    return ln_dos, energy
def U_F_S_ns(dbase, i):
    ln_dos, energy = dos_energy_ns(dbase, i)
    T_bins, dT, T_range = 1e3, 1/T_bins, np.arange(dT, 1, dT) # whats the max T?

    Z, U, F, S = ({} for k in range(4))
    

def U_F_S(dname):
    ln_dos, energy = dos_energy(dname + '-dos.dat')
    T_bins, dT, T_range = 1e3, 1.0/T_bins, np.arage(dT, 1.0, dT) 
    Quantities = {'U' = [], 'F' = [], 'S' = []}
    
    
    
            
    
