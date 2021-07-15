from __future__ import division, print_function
import numpy as np
from scipy.special import erfinv
import find_xi

def find_alpha(kT, sigma = 1, eps = 1):
    print(f'kT = {kT}')
    beta = 1/kT
    beta_eps = beta*eps
    p = (beta_eps**2*576 - beta_eps*3264 + 676)/147*sigma**12
    q = (beta_eps**3*27648 + beta_eps**2*55296 + beta_eps*67104 - 35152)/9261*sigma**18
    print(f'old p = {p}')
    print(f'old q = {q}')

    p = (-3*96/7*beta_eps - ((24*beta_eps-26)/7)**2)*sigma**12/3
    q = (2*((24*beta_eps-26)/7)**3 + 9*((24*beta_eps-26)/7)*96/7*beta_eps + 27*96/7*beta_eps)*sigma**18/27
    print(f'p = {p}')
    print(f'q = {q}')
    discriminant = -(4*p**3 + 27*q**2)
    print('discriminant = ', discriminant)

    if discriminant > 0:
        k = 1 # 0, 1, 2 are the three solutions
        print('thing in arccos', 3*q/(2*p)*np.sqrt(-3/p))
        t = 2*np.sqrt(-p/3)*np.cos(1/3*np.arccos(3*q/(2*p)*np.sqrt(-3/p)) - 2*np.pi*k/3)
        print(f't is {t}')
    elif p < 0:
        t = -2*abs(q)/q *np.sqrt(-p/3)*np.cosh(1/3*np.arccosh(-3*abs(q)/(2*p)*np.sqrt(-3/p)))
    else:
        t = -2*np.sqrt(p/3)*np.sinh(1/3*np.arcsinh(3*q/(2*p)*np.sqrt(3/p)))
    print('thing in sixth root', t - (beta_eps*24 - 26)*sigma**6/21)
    return (t - (beta_eps*24 - 26)*sigma**6/21)**(1/6)

def parameters(kT, sigma = 1, eps = 1):
    # alpha = ( 2.0/(1.0 + np.sqrt(np.log(2)*kT/eps)) )**(1/6)
    #alpha = find_alpha(kT, sigma, eps)
    alpha = find_xi.find_alpha(kT)
    # alpha = find_alpha(kT, sigma, eps)
    print(f'alpha = {alpha}')
    if np.isnan(alpha):
        exit(1)
    assert(alpha > 0)
    #Xi = alpha/(6*np.sqrt(np.pi))/(np.log(2) + np.sqrt(np.log(2)*eps/kT))
    Xi = find_xi.Xi(alpha, kT)

    diameter = 2**(1.0/6)*sigma
    return alpha, Xi, diameter

def V(r, sigma = 1, eps = 1):
    return 4*eps*((sigma/r)**12 - (sigma/r)**6) + eps

def Vprime(r, sigma = 1, eps = 1):
    return -4*eps*(12*(sigma/r)**12 - 6*(sigma/r)**6)/r
