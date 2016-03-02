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

def F_hardsphere(dbase, N):
    fbase = dbase + '/absolute/'
    F = 0
    # loop over files in the ./absolute/ directory
    for j in range(0, 400):
        valid = 0
        total = 0
        ratio = 0
        fname = fbase + '%05d' % (j)
        # if file exists, load the text from the file
        if os.path.isfile(fname+'.dat'):
            #open file and read in total valid and failed checks of small cell
            with (fname+'.dat') open as file:
                for line in file:
                    if 'valid small checks:' in line:
                        valid = line.split()[-1]
                    if 'total small checks:' in line:
                        total = line.split()[-1]
            # compute absolute F using read in data.
            ratio = valid / total
            F += -np.log(ratio)/ N
    return F

def U_F_S_ns(dbase, i):
    # Main function; take input of data directory and return all possible quantities
    ln_dos, energy = dos_energy_ns(dbase, i)
    Ns = eval(ln_dos.keys())
    T_bins, dT, T_range = 1e3, 1/T_bins, np.arange(dT, 1, dT) # change Tmax # Generate array of temperatures
    Z = np.zeros(len(T_range), len(Ns))
    Zinf, U, F, S = (np.zeros_like(Z) for i in range(4))
    
    for N in set(Ns): # Set iteration readability is great, but defining separate index seems sloppy...
        j=0
        F_HS = F_hardsphere(dbase, N)
        for T in set(T_range):
            # all computed quantities are excess, with the exception of final entropy (includes ideal gas)
            k=0
            ln_dos_boltz = ln_dos[N] - energy[N]/T
            dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max()) #overflow/underflow issues, need to keep lndos reasonable
            Z[k,j] = sum(dos_boltz)  #indexed by [T, N]
            U[k,j] = sum[energy[N]*dos_boltz]/Z[k,j]
            F[k,j] = -T*np.log(Z[k,j])

            # make absolute by equating entropy at T=inf. note that F_{ex,HS} \prop  S_{ex,HS}
            # need to take ln of sum of e^{ln_dos} to find number of states
            
            S_SW = np.log(sum(np.exp(ln_dos[N]))) + 3/2*N # U(T) = 3/2NkT probably
            S_HS = F_HS
            # F_{ex,HS} = -TS_HS; F_SW = U - TS_{ex,SW}
            # lim_{T \rightarrow \infty} U_{ex, SW} = constant; \therefore F_{ex,SW}(T=\infty) = -TS_{ex, SW }
            
            S_SW -= (S_SW - S_HS)
            F_SW = -T*(S_SW + 3/2*N)
            # Now, F_{ex, SW} is the absolute free energy required for the configuration. So, make calculated F absolute!
            F[k,j] -= (F[k,j] - F_SW)

            k += 1
        j += 1
    return U, F, S
