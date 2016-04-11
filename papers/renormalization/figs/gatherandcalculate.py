#!/usr/bin/python2
from __future__ import division
import numpy as np
import math
import string
import os

def lndos_energy_Ns(dbase):
    ln_dos = {}
    energy = {}
    ln_dos_hist = [] # need to initialize this correctly; how many elements?
    Ns = []
    for N in range(2,100): #maybe require N as an argument?
        fname = '%s/N%03d/lv-data-dos.dat' % (dbase, N)
        if os.path.isfile(fname):
            ln_dos_hist = np.loadtxt(fname, ndmin=2)
            ln_dos[N] = ln_dos_hist[:,1]
            energy[N] = ln_dos_hist[:,0] #slicing may be wrong, need to look at this harder. Also, if a file doesn't exist, index is no longer =N...
            Ns.append(N)
        else:
            pass
    return ln_dos, energy, np.array(Ns)

def Sexc_hardsphere_Ns(dbase):
    Ns = []
    S = []
    for N in range(2,100): # maybe go higher?
        fbase = '%s/N%03d/absolute/' % (dbase, N)
        if os.path.isfile(fbase+'00000.dat'):
            try:
                thisS = Sexc_hardsphere(dbase, N)
                S.append(thisS)
                Ns.append(N)
            except:
                print 'no data for N =', N
        else:
            print 'no such file as', fbase+'00000.dat'
    return np.array(S), np.array(Ns)

def Sexc_hardsphere(dbase, N):
    fbase = '%s/N%03d/absolute/' % (dbase, N)
    S = 0
    # loop over files in the ./absolute/ directory
    j = 0
    while True:
        valid = 0
        total = 0
        ratio = 0
        fname = fbase + '%05d.dat' % (j)
        # if file exists, load the text from the file
        if os.path.isfile(fname):
            #open file and read in total valid and failed checks of small cell
            with open(fname) as file:
                for line in file:
                    if 'valid small checks:' in line:
                        valid = int(line.split()[-1])
                    if 'total checks of small cell:' in line:
                        total = int(line.split()[-1])
            # compute absolute S using read in data.
            ratio = valid / total
            S += -np.log(ratio)/ N
        else:
            # we have found the last data file and can stop now
            break
        j += 1
    return S

def Uexc(ln_dos, energy, Ts):
    Ns = ln_dos.keys()
    U = np.zeros((len(Ts), len(Ns)))
    for j in range(len(Ns)):
        N = Ns[j]
        for k in range(len(Ts)):
            T = Ts[k]
            ln_dos_boltz = ln_dos[N] - energy[N]/T
            # Subtract of ln_dos_boltz.max() to keep dos_boltz reasonable
            dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max())
            U[k,j] = sum(energy[N]*dos_boltz)/sum(dos_boltz)
    return U

def Fexc(dbase, ln_dos, energy, Ts):
    Ns = ln_dos.keys()
    F = np.zeros((len(Ts), len(Ns)))
    for j in range(len(Ns)):
        N = Ns[j]
        Sexc_HS = Sexc_hardsphere(dbase, N)
        for k in range(len(Ts)):
            T = Ts[k]
            ln_dos_boltz = ln_dos[N] - energy[N]/T
            # Subtract of ln_dos_boltz.max() to keep dos_boltz reasonable
            dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max())
            Z = sum(dos_boltz)
            Zinf = sum(np.exp(ln_dos[N] - ln_dos_boltz.max()))
            F[k,j] = -T*Sexc_HS - T*np.log(Z/Zinf)
    return F

def Sexc(Uexc, Fexc, Ts):
    S = np.zeros_like(Uexc)
    for k in range(len(Ts)):
        T = Ts[k]
        S[k,:] = (Uexc[k,:] - Fexc[k,:])/T
    return S

# def U_F_S(dbase, i):
#     # Main function; take input of data directory and return all possible quantities
#     ln_dos, energy = dos_energy(dbase, i)
#     Ns = eval(ln_dos.keys())
#     T_bins, dT, T_range = 1e3, 1/T_bins, np.arange(dT, 1, dT) # change Tmax # Generate array of temperatures
#     Z = np.zeros(len(T_range), len(Ns))
#     Zinf, U, F, S = (np.zeros_like(Z) for i in range(4))
    
#     for j in range(len(Ns)): # Set iteration readability is great, but defining separate index seems sloppy...
#         N = Ns[j]
#         F_HS = F_hardsphere(dbase, N)
#         for k in range(len(T_range)):
#             # all computed quantities are excess, with the exception of final entropy (includes ideal gas)
#             T = T_range[k]
#             ln_dos_boltz = ln_dos[N] - energy[N]/T
#             dos_boltz = np.exp(ln_dos_boltz - ln_dos_boltz.max()) #overflow/underflow issues, need to keep lndos reasonable
#             Z[k,j] = sum(dos_boltz)  #indexed by [T, N]
#             U[k,j] = sum(energy[N]*dos_boltz)/Z[k,j]
#             F[k,j] = -T*np.log(Z[k,j]) fixme need max here? also ideal gas?

#             # make absolute by equating entropy at T=inf. note that F_{ex,HS} \prop  S_{ex,HS}
#             # need to take ln of sum of e^{ln_dos} to find number of states
            
#             S_SW = np.log(sum(np.exp(ln_dos[N]))) fixme deal with max of ln_dos
#             S_HS = F_HS
#             # F_{ex,HS} = -TS_HS; F_SW = U - TS_{ex,SW}
#             # lim_{T \rightarrow \infty} U_{ex, SW} = constant; \therefore F_{ex,SW}(T=\infty) = -TS_{ex, SW }
            
#             S_SW -= (S_SW - S_HS)
#             F_SW = -T*(S_SW + 3/2*N)
#             # Now, F_{ex, SW} is the absolute free energy required for the configuration. So, make calculated F absolute!
#             F[k,j] -= (F[k,j] - F_SW)

#     return U, F, S
