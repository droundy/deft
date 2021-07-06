from __future__ import division, print_function
#This program finds xi for a given temperture T (really kT) by setting the second virial coefficient
#computed with V_wca equal to the second virial coefficient computed with V_erf.
#It is referenced by wca_erf.py which is referenced by potential-plot.py and by w2-comparison.py.

from scipy import special
from scipy.special import erf
import numpy as np
import math

#Old Xi derived from derivatives (referenced here if want to revert to old Xi for comparison)------
#def Xi(alpha, kT, sigma = 1, eps = 1):
#   return alpha/(6*np.sqrt(np.pi))/(np.log(2) + np.sqrt(np.log(2)*eps/kT))
#--------------------------------------------------------------------------------------------------

#Compute B2_erf analytically (default method): 
def B2_erf_analytical(alpha, Xi, T) :
    return np.pi/3*((alpha**3 + 1.5*alpha*Xi**2)*(1+erf(alpha/Xi)) + 1/np.sqrt(np.pi)*(alpha**2*Xi + Xi**3)*np.exp(-(alpha/Xi)**2))


#Compute B2_wca numerically by evaluating the integral (default method):
def Vwca(r):
    return 4*((1/r)**12-(1/r)**6) +  1
    
def f_wca(r, T):
    return np.exp(-Vwca(r)/T) - 1

def B2_wca_numerical(T):
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    f = f_wca(r, T)
    f[0] = 1
    return -0.5*(4*np.pi*r**2*r[1]*f).sum()  

#Find Xi at a given  T (that is, kT)
def Xi(alpha, T, sigma = 1, eps = 1):
    B2wca = B2_wca_numerical(T)   #use with B2 wca numerical method  (default)
    xi_lo = 0
    xi_hi = sigma
    while xi_hi - xi_lo > 0.000001:
        xi_mid = 0.5*(xi_hi + xi_lo)
        if B2_erf_analytical(alpha, xi_mid, T) > B2wca:    #use with B2 erf analytical method (default)
            xi_hi = xi_mid
        else:
            xi_lo = xi_mid
    print('Xi is', xi_mid)
    return xi_mid

# # #CHECK---------------------------------------------------------------------
# #Find Xi at a given T   
# T=0.2575
# alpha=(2.0/(1+np.sqrt(T*np.log(2))))**(1.0/6)
# B2wca = B2_wca_numerical(T)   #use with B2 wca numerical method  (default)
# xi_lo = 0
# xi_hi = 1
# while xi_hi - xi_lo > 0.000001:
    # xi_mid = 0.5*(xi_hi + xi_lo)
    # if B2_erf_analytical(alpha, xi_mid, T) > B2wca:    #use with B2 erf analytical method (default)
        # xi_hi = xi_mid
    # else:
        # xi_lo = xi_mid
# Xi_at_T=xi_mid
# print(Xi_at_T, T)

