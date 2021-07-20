from __future__ import division, print_function
#This program finds xi for a given temperture T (really kT) by setting the second virial coefficient
#computed with V_wca equal to the second virial coefficient computed with V_erf.
#It is referenced by wca_erf.py which is referenced by potential-plot.py and by w2-comparison.py.

from scipy import special
from scipy.special import erf
import numpy as np
import math

#Compute B2_erf analytically (default method): 
def B2_erf_analytical(alpha, Xi, T) :
    return np.pi/3*((alpha**3 + 1.5*alpha*Xi**2)*(1+erf(alpha/Xi)) + 1/np.sqrt(np.pi)*(alpha**2*Xi + Xi**3)*np.exp(-(alpha/Xi)**2))


#Compute B2_wca numerically by evaluating the integral (default method):
def Vwca(r):
    return 4*((1/r)**12-(1/r)**6) +  1
    
def Vwca_prime(r):
	return 4*(-12*(1/r)**13+6*(1/r)**7)
    
def f_wca(r, T):
    return np.exp(-Vwca(r)/T) - 1
    
def fp_wca(r, T):
	return np.exp(-Vwca(r)/T)*(-Vwca_prime(r)/T)
	
def mean_fprime_wca_radius(T):    #mean of the values of r
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    fp = fp_wca(r, T)
    fp[0] = 0
    dr = r[1] - r[0]
    return (4*np.pi*r**3*dr*fp).sum()/(4*np.pi*r**2*dr*fp).sum() 
    
def mean_fprime_wca_radius_squared(T):
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    fp = fp_wca(r, T)
    fp[0] = 0
    dr = r[1] - r[0]
    return (4*np.pi*r**4*dr*fp).sum()/(4*np.pi*r**2*dr*fp).sum() 

def variance_fprime_wca_radius(T):
	return mean_fprime_wca_radius_squared(T) - (mean_fprime_wca_radius(T))**2
	
def find_Xi(alpha, T, sigma = 1, eps = 1):
	print("variance", variance_fprime_wca_radius(T))
	return (2*variance_fprime_wca_radius(T))**0.5

def B2_wca_numerical(T):
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    f = f_wca(r, T)
    f[0] = 1
    return -0.5*(4*np.pi*r**2*r[1]*f).sum()  
    
def find_alpha(T, sigma = 1, eps = 1):
    B2wca = B2_wca_numerical(T)   #use with B2 wca numerical method  (default)
    xi = find_Xi(0, T, sigma, eps)
    alpha_lo = 0
    alpha_hi = 3*sigma
    while alpha_hi - alpha_lo > 0.000001:
        alpha_mid = 0.5*(alpha_hi + alpha_lo)
        if B2_erf_analytical(alpha_mid, xi, T) > B2wca:    #use with B2 erf analytical method (default)
            alpha_hi = alpha_mid
        else:
            alpha_lo = alpha_mid
    print('Alpha is', alpha_mid, 'B2wca is', B2wca)
    return alpha_mid
    

