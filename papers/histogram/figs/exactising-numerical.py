#!/usr/bin/env python

from __future__ import division
import sys, os
import numpy as np
import math as mh

import matplotlib.pyplot as plt
import matplotlib

#------------------------------------------------------#
#------- COMPUTE THE EXACT REFERENCE FROM BEALE -------#

#LINK --> https://spot.colorado.edu/~beale/IsingExactMathematica.html

# define a Temperature Range
dT = 0.01
T = np.arange(1,100,dT)

# define parameter 'x' in terms of the Temperature
x = np.exp(-2/T)

# define Beta (not to be confused with 1/kBT) as per Paul D. Beale
b = 2*x - 2*x*x*x

def CalculateCoefSum(n,m,k):
    # define a[k_] as per Paul D. Beale
    a_coef = (1 + x*x)**2 - b*np.cos(np.pi*k/n)
    c2_coef_sum = 0
    for j in range(0, m+2, 2): # should include endpoint so m --> m+2
        c2_coef_sum += (mh.factorial(m)/mh.factorial(j)/mh.factorial(m-j))*(a_coef**2 - b*b)**(j//2) * a_coef**(m - j)
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
            c2_coef *= 2**(1 - 2*m)*((c2_coef_sum) + b**m)
            s2_coef *= 2**(1 - 2*m)*((c2_coef_sum) - b**m)

        z1 = 2**(m*n//2 - 1)*c2_coef
        z2 = 2**(m*n//2 - 1)*s2_coef

        c2_coef = 1
        s2_coef = 1
        for k in range(2, n , 2): # should include endpoint so n - 2 --> n
            c2_coef_sum = CalculateCoefSum(n,m,k)
            # define c2[k_] and s2[k_] as per Paul D. Beale
            c2_coef *= 2**(1 - 2*m)*((c2_coef_sum) + b**m)
            s2_coef *= 2**(1 - 2*m)*((c2_coef_sum) - b**m)

        z3 = 2**(m*n//2 - 1)*c0*cn*c2_coef
        z4 = 2**(m*n//2 - 1)*s0*sn*s2_coef
    else:
        print('Cannot determine the Parition Function for (ODD) spin systems')
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

# The general expression for the entropy is given (via Helmholtz Energy):
# F := E - TS = - kB T ln(Z) --> S(E) = kBln(Z(Beta)) + E/T
# E = - d/d(Beta) ln(Z) or equivalently E = kB T^2 d/dT ln(Z)

def EvaluatePartitionFunction(n,m,Z,Temp,dT):
    # convert Z to lnZ and take a derivative
    # K = J/kBT
    lnZ = 2*n*m/Temp + np.log(Z.astype(float))
    #Z = np.exp(2*n*m/Temp)*Z

    E = Temp*Temp*np.gradient(lnZ,dT)
    S =  lnZ + E/Temp

    return S,E

n,m = 32,32

Z_long = CalculateRef(n,m)

S,E = EvaluatePartitionFunction(n,m,Z_long,T,dT)
plt.plot(E,S)
plt.xlabel('Energy')
plt.ylabel('ln(S(E))')

plt.show()

# (2) Output in same file format as yaml-reference .py that way I can compare to
#     exact solution without having to run this code often and I won't have to change
#     the comparison file.

def OutputFile():
    np.savetxt(dirname,
          np.c_[energy, lndos, ps],
          fmt = ('%.16g'),
          delimiter = '\t',
          #newline = '# seed:',
          #newline = '# well_width: %g' % (well_width),
          #newline = '# ff: ',
          #newline = '# N: ',
          #newline = '# walls: ',
          #newline = '# cell dimensions: (%g, %g, %g)' % (x,y,z),
          #newline = '# translation_scale: %g' % (translation_scale),
          #newline = '# energy_levels: ',
          #newline = '# min_T: %g' % (min_T),
          #newline = '# max_entropy_state: ',
          #newline = '# min_important_energy: ',
          #newline = '# too_high_energy: %i' % (too_hi),
          #newline = '# too_low_energy: %i' % (too_lo),
          #newline = '',
          ##newline = '# WL Factor: %g' % (gamma),
          ##newline = '# iterations: %i' % (),
          ##newline = '# working moves: %i' % (),
          ##newline = '# total moves: %i' % (),
          ##newline = '# acceptance rate: %i' % (),
          #newline = '',
          ##newline = '# converged state: %i' % (),
          #newline = '# converged temperature: ',
          #newline = '# energy\t lndos\t ps\t lndos_tm: ',
          #newline = '\n version: created with yaml\n',
          header = 'comparison reference file\t(generated with python %s \n max_entropy_state: %i \n min_important_energy: %i \n energy\t lndos\t\t ps\t ' % (' '.join(sys.argv),max_entropy_state,min_important_energy))
