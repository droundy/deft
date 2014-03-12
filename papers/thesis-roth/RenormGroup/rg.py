from __future__ import division
import numpy as np
# import Hughes
import SW
# from scipy import integrate
import integrate

# Author: Dan Roth
# Email: Daniel.Edward.Roth@gmail.com
# Date: January 2014

# The theory here comes from:
## 'Application of a renormalization-group threatment to SAFT-VR'
## E. Forte, F. Llovell, L. F. Vegas, J. P. Trusler, A. Galindo
## Journal of Chemical Physics 134, 154102 (2011)

# Important constants
# k_B = 3.16681539628059e-6 # Boltzmann's constant in Hartree/Kelvin
# angstrom = 1.8897261 # An angstrom in bohr
k_B = 1 # define Boltzman constant as unity, since we will be working in SW units

# Square Well Parameters: set them all to unity
# sigma = 3.0342*angstrom # hard-sphere diameter
sigma = 2 # HS diameter
epsilon = 1 # depth of well
lambdaSW = 1.5 # range of well
R = sigma/2 # HS radius

# Temp
# T = 300 # Kelvin
T = 0.1 # SW Temp

# nL is the density to which you integrate when integrating over fluctuations
# See Forte 2011, text after eqn 7
# set it at unity
nL = 1

# Averaging volume constants
L = (8.509*(lambdaSW - 1)**2 - 4.078*(lambdaSW - 1) + 4.914)*sigma # reference wavelength; eqn (56), Forte 2011

# Fluctuation wavelength
def lambda_D(i):
    return 2**i*L

# Averaging volume
def VD(i):
    return (lambda_D(i)/2)**3


########### Various methods of recursive functions ######################
# # set number of evaluations counter
# num_f_evaluations = 0

# # Free energy density
# # eqn (5) from Forte 2011
# def f(T,n,i):
#     global num_f_evaluations
#     # f0 is just Huges
#     result = Hughes.f(T,nL)

#     # increase counter
#     num_f_evaluations += 1

#     for k in range(i+1):
#         result += -k_B*T*np.log(ID(T,nL,k)/ID_ref(T,nL,k))/VD(k)
# #        result += dfi(T,nL,k)
# #        print 'i:',i
# #        print 'f:',result
#     return result

# Free energy density
# eqn (5) from Forte 2011
# def frecursive(T,n,i):
#     if i == 0:
#         return Hughes.f(T,nL)
#     fiminusone = frecursive(T,n,i-1)
#     dfi = -k_B*T*np.log(ID(T,nL,i)/ID_ref(T,nL,i))/VD(i) # eqn (7), Forte 2011
#     return fiminusone + dfi

# Another option for defining stuff
# def f1(T, n):
#     return frecursive(T,nL,1)

# # incremental change in free energy density
# # eqn (7), Forte 2011
# def dfi(T,n,i):
#     return -k_B*T*np.log(ID(T,nL,i)/ID_ref(T,nL,i))/VD(i)

# Free energy density
def frecursive(T,n,i):
    # eqn (55) Forte 2011:
#    fnaught = Hughes.fid(T,n) + Hughes.fhs(T,n) + Hughes.fassoc(T,n) + Hughes.a2(n) # Hughes' a2 is the same as Forte's f2
    fnaught = SW.fid(T,n) + fhs(T,n) + SW.a2(n) # SW (Hughes) a2 is the same as Forte's f2
    # eqn (5) from Forte 2011:
    for j in range(i+1): # The function range(x) only goes up to x-1; range(x+1) will include x
        if j == 0:
            f = fnaught
            print '    j0',j
            print '    f =',f
        else:
            dfi = -k_B*T*np.log(ID(integrand_ID,T,nL,j)/ID_ref(integrand_ID_refT,nL,j))/VD(j) # eqn (7), Forte 2011
            f += dfi
            print '    j not 0',j
            print '    f =',f
    return f
##########################################################################

# Hard-sphere free energy
# Copied from Hughes --> Hughes uses atomic units, I need SW units
def fhs(T,n):
    # Phi 1
    Phi1 = -n0(n)*np.log(1-n3(n))

    # Phi 2
    Phi2 = (n1(n)*n2(n))/(1-n3(n))

    # Phi 3
    Phi3 = n2(n)**3*((n3(n)+(1-n3(n))**2*np.log(1-n3(n)))/(36*np.pi*n3(n)**2*(1-n3(n))**2))

    # free energy
    return k_B*T*(Phi1 + Phi2 + Phi3)

# Fundamental measure densities
def n0(n):
    return n

def n1(n):
    return n*R

def n2(n):
    return n*R**2*4*np.pi

def n3(n):
    return n*R**3*(4/3)*np.pi

# Integral over amplitudes (x) of wave-packet of length lambda_D
# similar to sum over density fluctuations n_D
# eqn(8), Forte 2011
def integrand_ID(x):
    return np.exp(-VD(i)/k_B/T*(fbarD(T,nL,x,i) + ubarD(T,nL,x,i)))
def ID(integrand,T,n,i):
    value,error = integrate.trapezoid(integrand,0,nL,1000)
    return value

# Reference for ID
# Evaluated at small enough wavelengths that UbarD should be negligible
# eqn (9), Forte 2011
def integrand_ID_ref(x):
    return np.exp(-VD(i)/k_B/T*fbarD(T,nL,x,i))
def ID_ref(integrand,T,n,i):
    value,error = integrate.trapezoid(integrand,0,nL,1000)
    return value

# Average within the considered subdomain, based on x
# eqn (10), Forte 2011
# For fluctuations in [0,1] (rather than [0,nL]), do nL(1 +/- x) instead of simply (nL +/- x), as it is in Forte
# This amounts to defining x* = x/nL --> nL-x = nL-(nL x*) = nL(1-x*); then simply writing x (instead of x*) in the program
def fbarD(T,n,x,i):
    value = (frecursive(T,nL*(1+x),i-1) + frecursive(T,nL*(1-x),i-1))/2 - frecursive(T,nL,i-1)
    return value

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,n,x,i):
    value = (u(T,nL*(1+x),lambda_D,i) + u(T,nL*(1-x),lambda_D,i))/2 - u(T,nL,lambda_D,i)
    return value

# Important values for u
# eqns (52-54), Forte 2011
alpha = 2/3*np.pi*epsilon*sigma**3*(lambdaSW**3 - 1)
omega2 = 1/5*sigma**2*(lambdaSW**5 - 1)/(lambdaSW**3 - 1)
gamma = 1/70*sigma**4*(lambdaSW**7 - 1)/(lambdaSW**3 - 1)

# Potential energy density for square well
# Eqn (51), Forte 2011
## The paper also includes m, which refers to the number of segments forming a chain; we do not deal with chains, so m = 1
def u(T,n,lambda_d,i):
    value = -(n)**2*(SW.gHS_eff(n)*alpha - (4*np.pi**2)/(2**(2*i+1)*L**2)*alpha*omega2 + ((4*np.pi**2)/(2**(2*i+1)*L**2))**2*alpha*gamma)
    return value

# Total
# Eqn (32), Forte 2011; f1 = fatt (see paragraph after eqn 55)
def ftot(T,n,i):
    return frecursive(T,n,i) + SW.a1SW(n) # SW.a1SW is the same as Forte's f1

# Grand free energy per volume
def phi(T,n,nparticular,i):
    mu = df_dn(T,nparticular,i)
    return ftot(T,n,i) - mu*n

# Derivative of free energy wrt number at const temp
def df_dn(T,n,i):
    # step size
    dn = 1e-6*n
    return (ftot(T,n + dn/2,i) - ftot(T,n - dn/2,i))/dn
