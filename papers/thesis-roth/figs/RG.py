from __future__ import division
import numpy as np
import SW
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
T = 0.1 # SW Temp

# nL is the density to which you integrate when integrating over fluctuations
# See Forte 2011, text after eqn 7
vol = 4/3*np.pi*R**3
nL = 0.3/vol

# Averaging volume constants
L = (8.509*(lambdaSW - 1)**2 - 4.078*(lambdaSW - 1) + 4.914)*sigma # reference wavelength; eqn (56), Forte 2011

# Fluctuation wavelength
def lambda_D(i):
    return 2**i*L

# Averaging volume
def VD(i):
    return (lambda_D(i)/2)**3

# Free energy density, defined by an iterative process
def fiterative(T,n,i):
    # eqn (55) Forte 2011:
    fnaught = SW.fid(T,n) + SW.fhs(T,n) + SW.a2(n) # SW (Hughes) a2 is the same as Forte's f2
    # eqn (5) from Forte 2011:
    for j in range(i+1): # The function range(x) only goes up to x-1; range(x+1) will include x
        if j == 0:
            f = fnaught
        else:
            IDvalue = ID(integrand_ID,T,nL,j)
            ID_refvalue = ID_ref(integrand_ID_ref,T,nL,j)
            dfi = -k_B*T*np.log(IDvalue/ID_refvalue)/VD(j) # eqn (7), Forte 2011
            f += dfi
    return f

# Integral over amplitudes (x) of wave-packet of length lambda_D
# similar to sum over density fluctuations n_D
# eqn(8), Forte 2011
def integrand_ID(T,x,i):
    argument = np.exp(-VD(i)/k_B/T*(fbarD(T,x,i) + ubarD(T,x,i)))
    fbarDvalue = fbarD(T,x,i)
    ubarDvalue =  ubarD(T,x,i)
    return argument

def ID(integrand,T,n,i):
    value = integrate.midpoint(lambda x: integrand_ID(T,x,i),0,1,1000)
    return value

# Reference for ID
# Evaluated at small enough wavelengths that UbarD should be negligible
# eqn (9), Forte 2011
def integrand_ID_ref(T,x,i):
    return np.exp(-VD(i)/k_B/T*fbarD(T,x,i))

def ID_ref(integrand,T,n,i):
    value = integrate.midpoint(lambda x: integrand_ID_ref(T,x,i),0,1,1000)
    return value

# Average within the considered subdomain, based on x
# eqn (10), Forte 2011
# For fluctuations in [0,1] (rather than [0,nL]), do nL(1 +/- x) instead of simply (nL +/- x), as it is in Forte
# This amounts to defining x* = x/nL --> nL-x = nL-(nL x*) = nL(1-x*); then simply writing x (instead of x*) in the program
def fbarD(T,x,i):
#    value = (fiterative(T,nL*(1+x),i-1) + fiterative(T,nL*(1-x),i-1))/2 - fiterative(T,nL,i-1)
    iplusx = fiterative(T,nL*(1+x),i-1)
    iminusx = fiterative(T,nL*(1-x),i-1)
    nochangex = fiterative(T,nL,i-1)
    value = (iplusx + iminusx)/2 - nochangex
    return value

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,x,i):
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
    return fiterative(T,n,i) + SW.a1SW(n) # SW.a1SW is the same as Forte's f1

# Grand free energy per volume
def phi(T,n,nparticular,i):
    mu = df_dn(T,nparticular,i)
    return ftot(T,n,i) - mu*n

# Derivative of free energy wrt number at const temp
def df_dn(T,n,i):
    # step size
    dn = 1e-6*n
    return (ftot(T,n + dn/2,i) - ftot(T,n - dn/2,i))/dn

# print 'ftot, i=1',ftot(0.1,0.1,1)
# print 'phi =',phi(0.1,0.1,0.1,0),'i=0'
# print 'phi =',phi(0.1,0.1,0.1,1),'i=1'
