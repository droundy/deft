#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, iv
from scipy.integrate import quad

def alpha(T):
    return (2/(1+np.sqrt(T*np.log(2))))**(1.0/6)

def B2_erf_numerical(Xi, T) :
    r = np.linspace(0, 10*alpha(T), 10000)
    f=(1.0/2)*(erf((r-alpha(T))/(Xi/np.sqrt(2)))-1)   #fixed sqrt(2)
    return -0.5*(4*np.pi*r**2*r[1]*f).sum()
B2_erf_numerical = np.vectorize(B2_erf_numerical)

def B2_erf_integrand(r, Xi, T) :
    erf_mayer_function=(1.0/2)*(erf((r-alpha(T))/(Xi/np.sqrt(2)))-1)   #fixed sqrt(2)
    return (-1.0/2)*(4*np.pi)*(erf_mayer_function)*r*r
def B2_erf_with_quad(Xi, T):
    val = quad(B2_erf_integrand, 0, np.inf, args=(Xi, T))
    # print('quad says', val)
    return val[0]

def b2(xi, T):   #Analytical method ERROR here...
    a = alpha(T)
    return np.pi/3*( (a**3 + 1.5*a*xi**2)*(1+erf(a/xi)) + 1/np.sqrt(np.pi)*(a**2*xi + xi**3)*np.exp(-(a/xi)**2) )    #fix sqrt(2) !

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

xi = np.linspace(0.0001, 1, 1000)

for T in [0.1, 0.5,1, 200]:
    plt.figure()
    plt.plot(xi, b2(xi, T), label=r'$B_2(\Xi)$ analytical')
    plt.plot(xi, B2_erf_numerical(xi, T), 'b--', label=r'$B_2(\Xi)$ numerical')
    with_quad = np.zeros_like(xi)
    for j in range(len(xi)):
        with_quad[j] = B2_erf_with_quad(xi[j], T)
    plt.plot(xi, with_quad, 'c:', label=r'$B_2(\Xi)$ quad')
    plt.axvline(alpha(T), color='xkcd:purple', linestyle=':', label=r'$\alpha$')
    plt.axhline(b2wca(T), color='g', label=r'$B_2(WCA)$')
    # plt.axhline(b2wca_bessel(T), color='r')
    plt.xlabel(r'$\Xi$')
    plt.ylabel(r'$B_2$')
    plt.title('T={}'.format(T))
    plt.legend(loc='best')


    # r = np.linspace(0.01, rmax, 1000)
    # plt.figure()
    # plt.plot(r, f_wca(r, T))
    #plt.plot(r, Vwca(r))

# T = np.arange(0.5, 2.0, 0.1)
# plt.figure()
# b2 = np.zeros_like(T)
# b2b = np.zeros_like(T)
# for i in range(len(T)):
#     b2[i] = b2wca(T[i])
#     b2b[i] = b2wca_bessel(T[i])
# plt.plot(T, b2)
# plt.plot(T, b2b)
plt.show()
