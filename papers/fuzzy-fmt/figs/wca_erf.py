from __future__ import division, print_function
import numpy as np
from scipy.special import erfinv
import find_xi

def find_alpha(kT, sigma = 1, eps = 1):
    c = 4*eps/kT
    B = 1/(24*c)
    p = (3*7*B-1)/3
    q = 1/27*(-2+9*B*(6*c-27) + 27*7*B)
    k = 2
    print(f'p = {p}')
    print(f'q = {q}')
    return 2*np.sqrt(-p/3)*np.cos((1/3)*np.arccos(3*q*np.sqrt(-3/p)/(2*p) + 2*np.pi*k/3)) - 1/3

def parameters(kT, sigma = 1, eps = 1):
    alpha = ( 2.0/(1.0 + np.sqrt(np.log(2)*kT/eps)) )**(1/6)
    # alpha = find_alpha(kT, sigma, eps)
    print(f'alpha = {alpha}')
    assert(alpha > 0)
    #Xi = alpha/(6*np.sqrt(np.pi))/(np.log(2) + np.sqrt(np.log(2)*eps/kT))
    Xi = find_xi.Xi(alpha, kT)

    diameter = 2**(1.0/6)*sigma
    return alpha, Xi, diameter

def V(r, sigma = 1, eps = 1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6) + eps

def Vprime(r, sigma = 1, eps = 1):
    return -4*eps*(12*(sigma/r)**12 - 6*(sigma/r)**6)/r
