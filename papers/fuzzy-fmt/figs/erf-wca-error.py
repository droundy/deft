from __future__ import division
from numpy import *
import matplotlib.pyplot as plt
from scipy.special import erf

min_temp = 0.1
max_temp = 1.3
eps = 1
R = 1
n = 0.5
sigma = R*2**(5.0/6.0)
Temps = arange(min_temp, max_temp, 0.1)
alphas = sigma*(2/(1 + sqrt(Temps *log(2)/eps)))**(1.0/6.0)
Xis = alphas/(6*sqrt(pi)*(log(2) + sqrt(Temps *log(2)/eps)))

N = 1000

def V_wca(r):
    if r < 2*R:
        return 4*eps*((sigma/r)**12 - (sigma/r)**6) + eps
    else:
        return 0

def V_erf(r,kT, a, X):
    return -kT*log((erf((r-a)/X)+1)/2)

plt.figure(1)
for i in range(len(Temps)):
    x = linspace(2*R, alphas[i]*1.5,N)
    dx = (1.3*alphas[i]-2*R)/N
    error = zeros_like(x)

    for j in arange(1,len(x)):
        error[j] = error[j-1] + 4*pi*n*( V_erf(x[j],Temps[i], alphas[i], Xis[i]))*x[i]**2*dx
    plt.plot(x, error, label = r'$kT/\epsilon = %g$' % Temps[i] )

plt.xlim(2*R,2.4)
plt.ylabel("$\delta E/\epsilon$")
plt.xlabel('$r/R$')
plt.title('Integrated difference $V_{wca}$ and $V_{erf}$ for r > 2R')
plt.legend(loc='best')

plt.figure(2)
for i in range(len(Temps)):
    x = linspace(alphas[i], alphas[i]*1.3,N)
    dx = alphas[i]*(1.3-1)/N
    error = zeros_like(x)

    for j in arange(1,len(x)):
        error[j] = error[j-1] + 4*pi*n**2*(V_wca(x[j]) -  V_erf(x[j],Temps[i], alphas[i], Xis[i]))*x[i]**2*dx
    plt.plot(x, error, label = r'$kT/\epsilon = %g$' % Temps[i] )

plt.xlim(alphas[-1],alphas[0]*1.15)
plt.ylabel("$\delta E/\epsilon$")
plt.xlabel('$r/R$')
plt.title('Integrated difference between $V_{wca}$ and $V_{erf}$')
plt.legend(loc='best')
plt.show()
