#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import sympy as sy
import mpmath as mp
import math as mh
from sympy import poly

import matplotlib.pyplot as plt
import matplotlib

#------------------------------------------------------#
#------- COMPUTE THE EXACT REFERENCE FROM BEALE -------#

#LINK --> https://spot.colorado.edu/~beale/IsingExactMathematica.html

x,T = sy.symbols('x T')

# global coefficient only depends on x
b = 2*x - 2*x*x*x

# define Floor function as implemented by Mathematica
#def Floor(num, multiple):
    #return num - (num % multiple)

#def RoundAccuracy(n,m):
    ## set the numerical precision according to paper
    ## (took forever to figure out what 'acc' was)!
    #acc = Floor(n*m*np.log(2)/np.log(10),1.5)
    #return acc

def CalculateCoefSum(n,m,k):
    # define a[k_] as per Paul D. Beale
    a_coef = (1 + x*x)**2 - b*mp.cos(np.pi*k/n)
    c2_coef_sum = 0
    for j in range(0, m+2, 2): # should include endpoint so m --> m+2
        c2_coef_sum += sy.expand((mh.factorial(m)/mh.factorial(j)/mh.factorial(m-j))*(a_coef**2 - b*b)**(j//2) * a_coef**(m - j))
    return c2_coef_sum

def CalculateRef(n, m):
    if (n % 2 == 0):
    
        # define coefficients (x = exp(-2K) = exp(-2J/kbT).
        # the coefficients c0,s0,cn,sn are just numbers contrary
        # to what the notation would have you believe.
        # dividing by 2**(m/2) is different than paper?
        c0 = ((1 - x)**m + (x*(1 + x))**m) / 2**(m//2)
        s0 = ((1 - x)**m - (x*(1 + x))**m) / 2**(m//2)
        cn = ((1 + x)**m + (x*(1 - x))**m) / 2**(m//2)
        sn = ((1 + x)**m - (x*(1 - x))**m) / 2**(m//2)
        
        c2_coef = 1
        s2_coef = 1
        for k in range(1, n , 2): # should include endpoint so n - 1 --> n
            c2_coef_sum = CalculateCoefSum(n,m,k)
            # define c2[k_] and s2[k_] as per Paul D. Beale
            c2_coef *= sy.expand(2**(1 - 2*m)*((c2_coef_sum) + b**m))
            s2_coef *= sy.expand(2**(1 - 2*m)*((c2_coef_sum) - b**m))
        
        z1 = 2**(m*n//2 - 1)*c2_coef
        z2 = 2**(m*n//2 - 1)*s2_coef
        
        c2_coef = 1
        s2_coef = 1
        for k in range(2, n , 2): # should include endpoint so n - 2 --> n
            c2_coef_sum = CalculateCoefSum(n,m,k)
            # define c2[k_] and s2[k_] as per Paul D. Beale
            c2_coef *= sy.expand(2**(1 - 2*m)*((c2_coef_sum) + b**m))
            s2_coef *= sy.expand(2**(1 - 2*m)*((c2_coef_sum) - b**m))
    
        z3 = 2**(m*n//2 - 1)*c0*cn*c2_coef
        z4 = 2**(m*n//2 - 1)*s0*sn*s2_coef
    else:
        print('Cannot determine Parition Function for (ODD) spin systems')
        ## check whether system is even or odd
        #if (n % 2 == 0):
            #z1 = 2^(m n/2 - 1) Product[c2[k], k, 1, n - 1, 2];
            #z2 = 2^(m n/2 - 1) Product[s2[k], {k, 1, n - 1, 2}];
            #z3 = 2^(m n/2 - 1) c[0] c[n] Product[c2[k], {k, 2, n - 2, 2}];
            #z4 = 2^(m n/2 - 1) s[0] s[n] Product[s2[k], {k, 2, n - 2, 2}], 
        #else: # we could include the cases of m != n
            ##z1 = 2^(m n/2 - 1) c[n] Product[c2[k], {k, 1, n - 1, 2}];
            ##z2 = 2^(m n/2 - 1) s[n] Product[s2[k], {k, 1, n - 1, 2}];
            ##z3 = 2^(m n/2 - 1) c[0] Product[c2[k], {k, 2, n - 1, 2}];
            ##z4 = 2^(m n/2 - 1) s[0] Product[s2[k], {k, 2, n - 1, 2}]]; 
    return z1 + z2 + z3 + z4


Z_long = sy.expand(CalculateRef(8,8))
Z_coeffs = poly(Z_long, x).all_coeffs()

Z = np.array(Z_coeffs,dtype=np.int64)

index = []
for i in range(len(Z_coeffs)):
    if (i % 2 != 0):
        index.append(i)

Z_new = np.delete(Z, index)
print(Z_new) # Test works for 8x8!

Temperature = np.arange(1,3,0.1)

def EvaluatePartitionFunction(n,m,Z,Temp):
    Z_ln_eval = []
    C_eval = []
    # convert Z and take a derivative
    Z = -sy.exp(-n*m*x)*Z
    Z = Z.subs(x,sy.exp(-2/T))
    Z =  sy.diff(1/T *Z,T)
    C = T*sy.diff(Z,T)
    for i in range(len(Temp)):
        Z_ln_eval.append(float(str(sy.re(sy.log(Z.subs(T,Temp[i]))))))
        C_eval.append(float(str(sy.re(C.subs(T,Temp[i])))))
    return Z_ln_eval,C_eval

S,C = EvaluatePartitionFunction(8,8,Z_long,Temperature)
#print('The entropy S is given:',S)
#print('The Heat Capacity C is given:',C)
plt.plot(1/Temperature,C)
plt.plot(1/Temperature,S)
plt.show()
## the Partition Function Z is then given by:
#Z = np.exp(2*n*m*K)*(z1 + z2 + z3 + z4)
##----------------------------------------------------------#

  

