from __future__ import division
import numpy as np

def parameters(kT, sigma = 1, eps = 1):
    alpha = ( 2.0/(1 + 6*np.sqrt(np.log(2)*kT)) )**(1.0/6)
    Xi = alpha/np.sqrt(np.pi)/(6*np.log(2) + np.sqrt(np.log(2)/kT))

    diameter = 2**(1.0/6)
    return alpha, Xi, diameter

def V(r, sigma = 1, eps = 1):
    return eps/9*((sigma/r)**12 - (sigma/r)**6) + eps/36

def Vprime(r, sigma = 1, eps = 1):
    return -eps/9*(12*(sigma/r)**12 - 6*(sigma/r)**6)/r
