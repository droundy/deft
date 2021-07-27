from __future__ import division, print_function
import numpy as np
from scipy.special import erfinv
import find_xi

def parameters(kT, sigma = 1, eps = 1):
    alpha = find_xi.find_alpha(kT)
    print(f'alpha = {alpha}')
    if np.isnan(alpha):
        exit(1)
    assert(alpha > 0)
    Xi = find_xi.find_Xi(alpha, kT)

    diameter = 2**(1.0/6)*sigma
    return alpha, Xi, diameter

def V(r, sigma = 1, eps = 1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6) + eps

def Vprime(r, sigma = 1, eps = 1):
    return -4*eps*(12*(sigma/r)**12 - 6*(sigma/r)**6)/r
