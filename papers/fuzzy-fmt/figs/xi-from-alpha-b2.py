#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, iv

def alpha(T):
    return (2/(1+np.sqrt(T*np.log(2))))**(1./6)

def b2(xi, T):
    a = alpha(T)
    return np.pi/3*( (a**3 + 1.5*a*xi**2)*(1+erf(a/xi)) + 1/np.sqrt(np.pi)*(a**2*xi + xi**3)*np.exp(-(a/xi)**2) )

def b2wca_bessel(T):
    return -2*np.sqrt(2)*np.pi**2/(6*T)*np.exp(-0.5/T)*(
        iv(-0.75, 0.5/T) + iv(0.75, 0.5/T) - iv(0.25, 0.5/T) - iv(-0.25, 0.5/T)
    )

def Vwca(r):
    return 4*((1/r)**12-(1/r)**6) +  1

rmax = 2.0**(1/6)

print(Vwca(rmax))

def f_wca(r, T):
    return np.exp(-Vwca(r)/T) - 1

def b2wca(T):
    r = np.linspace(0, rmax, 10000)
    f = f_wca(r, T)
    f[0] = 1
    return -0.5*(4*np.pi*r**2*r[1]*f).sum()

xi = np.linspace(0, 1, 1000)

for T in [0.1, 0.5,1, 200]:
    plt.figure()
    plt.plot(xi, b2(xi, T), label=r'$\Xi$')
    plt.axvline(alpha(T))
    plt.axhline(b2wca(T), color='g')
    plt.axhline(b2wca_bessel(T), color='r')
    plt.xlabel(r'$\Xi$')
    plt.ylabel(r'$B_2$')
    plt.title(T)


    # r = np.linspace(0.01, rmax, 1000)
    # plt.figure()
    # plt.plot(r, f_wca(r, T))
    #plt.plot(r, Vwca(r))

T = np.arange(0.5, 2.0, 0.1)
plt.figure()
b2 = np.zeros_like(T)
b2b = np.zeros_like(T)
for i in range(len(T)):
    b2[i] = b2wca(T[i])
    b2b[i] = b2wca_bessel(T[i])
plt.plot(T, b2)
plt.plot(T, b2b)
plt.show()
