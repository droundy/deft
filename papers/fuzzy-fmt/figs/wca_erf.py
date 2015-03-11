from __future__ import division
import numpy as np
from scipy.special import erfinv

def parameters(kT, sigma = 1, eps = 1):
    alpha = ( 2.0/(1.0 + np.sqrt(np.log(2)*kT/eps)) )**(1.0/6)
    Xi = alpha/(6*np.sqrt(np.pi))/(np.log(2) + np.sqrt(np.log(2)*eps/kT))

    diameter = 2**(1.0/6)*sigma
    return alpha, Xi, diameter

def V(r, sigma = 1, eps = 1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6) + eps

def Vprime(r, sigma = 1, eps = 1):
    return -4*eps*(12*(sigma/r)**12 - 6*(sigma/r)**6)/r
