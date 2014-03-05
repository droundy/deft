from __future__ import division
import numpy as np
import Hughes
from scipy import integrate

# Author: Dan Roth
# Email: Daniel.Edward.Roth@gmail.com
# Date: January 2014

# The theory here comes from:
## 'Application of a renormalization-group threatment to SAFT-VR'
## E. Forte, F. Llovell, L. F. Vegas, J. P. Trusler, A. Galindo
## Journal of Chemical Physics 134, 154102 (2011)

# Important constants
k_B = 3.16681539628059e-6 # Boltzmann's constant in Hartree/Kelvin
angstrom = 1.8897261 # An angstrom in bohr

# Square Well Parameters: set them all to unity
# sigma = 3.0342*angstrom # hard-sphere diameter
sigma = 1 # hard sphere diameter
epsilon = 1 # Depth of square well
lambdaSW = 1.5 # Length (range) of square well

# wave packet amplitude
x = 1 # ???

# Temp
T = 300 # Kelvin

# Density
nD = 1 # shortest wavelength fluctiations
nl = 1 # everything else
ns = nD + nl # total density

# Use this for now
nl = 4.9388942e-3 # atoms bohr^-3

# Averaging volume constants
i = 0 # iteration counter; start at naught
L = (8.509*(lambdaSW - 1)**2 - 4.078*(lambdaSW - 1) + 4.914)*sigma # reference wavelength; eqn (56), Forte 2011

# Fluctuation wavelength
def lambda_D(i):
#    print 'lambda_D:',2**i*L
    return 2**i*L

# Averaging volume
def VD(i):
#    print 'VD:',(lambda_D(i)/2)**3
    return (lambda_D(i)/2)**3

# set number of evaluations counter
num_f_evaluations = 0

# # Free energy density
# # eqn (5) from Forte 2011
# def f(T,n,i):
#     global num_f_evaluations
#     # f0 is just Huges
#     result = Hughes.f(T,nl)

#     # increase counter
#     num_f_evaluations += 1

#     for k in range(i+1):
#         result += -k_B*T*np.log(ID(T,nl,k)/ID_ref(T,nl,k))/VD(k)
# #        result += dfi(T,nl,k)
# #        print 'i:',i
# #        print 'f:',result
#     return result

# Free energy density
# eqn (5) from Forte 2011
# def frecursive(T,n,i):
#     if i == 0:
#         return Hughes.f(T,nl)
#     fiminusone = frecursive(T,n,i-1)
#     dfi = -k_B*T*np.log(ID(T,nl,i)/ID_ref(T,nl,i))/VD(i) # eqn (7), Forte 2011
#     return fiminusone + dfi

# Another thought
def frecursive(T,n,i):
    fnaught = Hughes.fid(T,n) + Hughes.fhs(T,n) + Hughes.fassoc(T,n) + Hughes.a2(n) # I'm pretty sure that Hughes' a2 is the same as Forte's f2
    for j in range(i):
        if j == 0:
            f = fnaught
        else:
            dfi = -k_B*T*np.log(ID(T,nl,i)/ID_ref(T,nl,i))/VD(i) # eqn (7), Forte 2011
            f += dfi
    return f

# Another option for defining stuff
# def f1(T, n):
#     return frecursive(T,nl,1)

# # incremental change in free energy density
# # eqn (7), Forte 2011
# def dfi(T,n,i):
#     return -k_B*T*np.log(ID(T,nl,i)/ID_ref(T,nl,i))/VD(i)

# Integral over amplitudes (x) of wave-packet of length lambda_D
# similar to sum over density fluctuations n_D
# eqn(8), Forte 2011
def ID(T,n,i):
    def integrand(x):
        return np.exp(-VD(i)/k_B/T*(fbarD(T,nl,x,i) + ubarD(T,nl,x,i)))
    value,error = integrate.quad(integrand,0,nl)
#    print 'ID: ',value
    return value

# Reference for ID
# Evaluated at small enough wavelengths that UbarD should be negligible
# eqn (9), Forte 2011
def ID_ref(T,n,i):
    def integrand(x):
        return np.exp(-VD(i)/k_B/T*fbarD(T,nl,x,i))
    value,error = integrate.quad(integrand,0,nl)
#    print 'ID_ref: ',value
    return value

# Average within the considered subdomain, based on x
# eqn (10), Forte 2011
def fbarD(T,n,x,i):
    value = (frecursive(T,nl+x,i-1) + frecursive(T,nl-x,i-1))/2 - frecursive(T,nl,i-1)
#    print 'fbarD:',value
    return value

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,n,x,i):
    value = (u(T,nl+x,lambda_D) + u(T,nl-x,lambda_D))/2 - u(T,nl,lambda_D)
#    print 'UbarD:',value
    return value

# Important values for u
# eqns (52-54), Forte 2011
alpha = 2/3*np.pi*epsilon*sigma**3*(lambdaSW**3 - 1)
omega2 = 1/5*sigma**2*(lambdaSW**5 - 1)/(lambdaSW**3 - 1)
gamma = 1/70*sigma**4*(lambdaSW**7 - 1)/(lambdaSW**3 - 1)

# Potential energy density for square well
# Eqn (51), Forte 2011
## The paper also includes m, which refers to the number of segments forming a chain; we do not deal with chains, so m = 1
def u(T,n,lambda_d):
    value = -(n)**2*(Hughes.gHS_eff(n_eff)*alpha - (4*np.pi**2)/(2**(2*i+1)*L**2)*alpha*omega2 + ((4*np.pi**2)/(2**(2*i+1)*L**2))**2*alpha*gamma)
#    print 'U:',value
    return value

# print f(T,nl,1)
n_eff = nl
m=1
#print frecursive(T,nl,1)

# Total
# Eqn (32), Forte 2011; f1 = fatt (see paragraph after eqn 55)
def ftot(T,n,i):
    return frecursive(T,n,i) + Hughes.a1SW(n) #f1(T,n)
