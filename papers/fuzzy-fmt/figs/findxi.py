import numpy as np
import scipy.special

def find_alpha(T) :
    return (2.0/(1+np.sqrt(T*np.log(2))))**(1.0/6)

def B2_erf_analytical(Xi,T) :
    alpha = find_alpha(T)
    return np.pi/3*((alpha**3 + 1.5*alpha*Xi**2)*(1+scipy.special.erf(alpha/Xi)) + 1/np.sqrt(np.pi)*(alpha**2*Xi + Xi**3)*np.exp(-(alpha/Xi)**2))

def Vwca(r):
    return 4*((1/r)**12-(1/r)**6) +  1

def B2_wca_numerical(T):
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    dr = r[1]
    f = np.exp(-Vwca(r)/T) - 1
    f[0] = 1
    return -0.5*(4*np.pi*r**2*dr*f).sum()

def find_Xi(T):
    B2wca = B2_wca_numerical(T)
    xi_lo = 0
    xi_hi = 1
    while xi_hi - xi_lo > 0.000001:
        xi_mid = 0.5*(xi_hi + xi_lo)
        if B2_erf_analytical(xi_mid, T) > B2wca:
            xi_hi = xi_mid
        else:
            xi_lo = xi_mid
    return xi_mid
