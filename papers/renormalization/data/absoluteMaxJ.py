#!/usr/bin/python

from __future__ import division
from math import pi       # REALLY don't need all of math
import os, numpy as np, sys

#assert(not os.system("fac ../../../free-energy-monte-carlo ../../../free-energy-monte-carlo-infinite-case"))
def free_energy_over_kT_to_ff(free_energy): # Find an ff that corresponds to given free energy
	    # bisection
	    n = 0
	    low,high = 0,1
	    resolution = 1E-8
	    while n < 500: #better value?
	        ff_guess = (low+high)/2
	        free_energy_next =  (4*ff_guess - 3*ff_guess**2)/(1-ff_guess)**2
	        if abs(free_energy_next - free_energy) < resolution:
	            return ff_guess
	        elif free_energy_next > free_energy:
	            high = (high+low)/2
	        else:
	            low = (high +low)/2
	        n+=1
	    return ff_guess
def findMaxJ(i,ww,L,N):
    R=1
    L_i = L*(2**i)
    #incorperation of RG parameter i
    #print("\n----------------------FREE ENERGY MONTE CARLO------------------------------\n\n")
    # The maximum free energy is near the point where
    # Carnahan-Starling predict a filling fraction of 0.7, which is
    # not reasonable to simulate.
    approximate_free_energies = np.arange(np.log(2), 3*N, np.log(2))/N
    ffs = np.zeros_like(approximate_free_energies)
    for k in xrange(len(ffs)):
        ffs[k] = free_energy_over_kT_to_ff(approximate_free_energies[k])
    ff_goal = (4*np.pi/3*R**3)*N/(L_i)**3 # this is the density
    maxJ=0
    for j in xrange(len(ffs)-2):
        maxJ+=1
        if ffs[j+2]>ff_goal: break
    return maxJ
