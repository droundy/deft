from __future__ import division
import numpy as np

# The following is not from the standard library. Find it in
# (...)/deft/papers/thesis-roth/figs
import integrate

#####################################################################
# Author: Dan Roth                                                  #
# Email: Daniel.Edward.Roth@gmail.com                               #
# Date: June 2014                                                   #
#                                                                   #
# This is the second-ish version. Changes:                          #
## Re-written for clarity                                           #
## Work directly with mu for coexistence                            #
## Better comments                                                  #
## Organized references at end of code                              #
#                                                                   #
# The theory comes mostly from Reference (1)                        #
#####################################################################

# Important constants
k_B = 1 # define Boltzman constant as unity, since we will be working in arbitrary units

# Square Well Parameters
sigma = 2 # HS diameter
R = sigma/2 # HS radius
epsilon = 1 # depth of well
lambdaSW = 1.5 # range of well

# Important values for u
# eqns (52-54), Forte 2011
alpha = 2/3*np.pi*epsilon*sigma**3*(lambdaSW**3 - 1)
omega2 = 1/5*sigma**2*(lambdaSW**5 - 1)/(lambdaSW**3 - 1)
gamma = 1/70*sigma**4*(lambdaSW**7 - 1)/(lambdaSW**3 - 1)

# Averaging volume
L = (8.509*(lambdaSW - 1)**2 - 4.078*(lambdaSW - 1) + 4.914)*sigma # reference wavelength; eqn (56), Forte 2011

# Fluctuation wavelength
def lambda_D(i):
    return 2**i*L

# Averaging volume, i.e. volume of Wilson Phase Cell
## See reference (2)
def VD(i):
    return (lambda_D(i)/2)**3

# Free energy density, defined by an iterative process
def fiterative(T,n,i):
    # eqn (55) Forte 2011:
    fnaught = SWfid(T,n) + SWfhs(T,n) + a2(n)/k_B/T*n # SW (and Hughes) a2/kT is the same as Forte's f2

    # eqn (5) from Forte 2011:
    for j in range(i+1): # The function range(y) only goes up to y-1; range(y+1) will include y
        if j == 0:
            f = fnaught
        else:
            IDvalue = ID(integrand_ID,T,n,j)
            ID_refvalue = ID_ref(integrand_ID_ref,T,n,j)
            dfi = -k_B*T*np.log(IDvalue/ID_refvalue)/VD(j) # eqn (7), Forte 2011
            f += dfi

    return f

# Integral over amplitudes (x) of wave-packet of length lambda_D
# similar to sum over density fluctuations n_D
# eqn(8), Forte 2011
def integrand_ID(T,n,x,i):
    argument = np.exp(-VD(i)/k_B/T*(fbarD(T,n,x,i) + ubarD(T,n,x,i)))
    return argument

def ID(integrand,T,n,i):
    maxn = 1/(sigma**3*np.pi/6)
    maxx = np.minimum(np.ones_like(n),maxn/n)
    return integrate.midpoint(lambda x: integrand_ID(T,n,x,i),0,maxx,100)

# Reference for ID
# Evaluated at small enough wavelengths that UbarD should be negligible
# eqn (9), Forte 2011
def integrand_ID_ref(T,n,x,i):
    return np.exp(-VD(i)/k_B/T*fbarD(T,n,x,i))

def ID_ref(integrand,T,n,i):
    maxn = 1/(sigma**3*np.pi/6)
    maxx = np.minimum(np.ones_like(n),maxn/n)
    return integrate.midpoint(lambda x: integrand_ID_ref(T,n,x,i),0,maxx,100)

# Average within the considered subdomain, based on x
## eqn (10), Forte 2011
## For fluctuations in [0,1] (rather than [0,n]), do n(1 +/- x) instead of simply (n +/- x), as it is in Forte
## This amounts to defining x* = x/n --> n-x = n-(n x*) = n(1-x*); then simply writing x (instead of x*) in the program
def fbarD(T,n,x,i):
    iplusx = fiterative(T,n*(1+x),i-1)
    iminusx = fiterative(T,n*(1-x),i-1)
    nochangex = fiterative(T,n,i-1)
    value = (iplusx + iminusx)/2 - nochangex
    return value

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,n,x,i):
    value = (u(T,n*(1+x),lambda_D,i) + u(T,n*(1-x),lambda_D,i))/2 - u(T,n,lambda_D,i)
    return value

# Potential energy density for square well
# Eqn (51), Forte 2011
## The paper also includes a factor m, which refers to the number of segments forming a chain
## we do not deal with chains, so m = 1
def u(T,n,lambda_d,i):
    value = -(n)**2*(gHS_eff(n)*alpha - (4*np.pi**2)/(2**(2*i+1)*L**2)*alpha*omega2 + ((4*np.pi**2)/(2**(2*i+1)*L**2))**2*alpha*gamma)
    return value

# Total free energy
# Eqn (32), Forte 2011; f1 = fatt (see paragraph after eqn 55)
def ftot(T,n,i):
    return fiterative(T,n,i) + a1SW(n)*n # a1SW*n is the same as Forte's f1

# Pressure
# See Jess Hughes' project report, equation (38)
def P(T,n,i):
  return n*df_dn(T,n) - ftot(T,n,i)

# Grand free energy per volume
def phi(T,n,mu,i):
    return ftot(T,n,i) - mu*n

def df_dn(T,n,i):
    dn = 1e-8*n
    return (ftot(T,n+dn/2,i) - ftot(T,n-dn/2,i))/dn

#######################################
# I define the following functions here, rather than calling from SW.py
# This so that you can use the functions without worring if (for example) HS radius is set differently between the two programs

# eta is the packing fraction; volume times number density
# In the homogeneous case, eta is the same as n3(n) (from FMT; see Hughes 2013)
def  eta(n):
    return n*R**3*(4/3)*np.pi

# effective packing fraction
# Gil-Villegas eqn (36)
def eta_eff(n):
    c1 = 2.25855-1.50349*lambdaSW+0.249434*lambdaSW**2
    c2 = -0.669270+1.40049*lambdaSW-0.827739*lambdaSW**2
    c3 = 10.1576-15.0427*lambdaSW+5.30827*lambdaSW**2
    return c1*eta(n) + c2*eta(n)**2 + c3*eta(n)**3

# derivative of eta_eff wrt eta
def deta_eff_deta(n):
    c1 = 2.25855-1.50349*lambdaSW+0.249434*lambdaSW**2
    c2 = -0.669270+1.40049*lambdaSW-0.827739*lambdaSW**2
    c3 = 10.1576-15.0427*lambdaSW+5.30827*lambdaSW**2
    return c1+2*c2*eta(n)+3*c3*eta(n)**2

# 1st term in High-temp perturbation expansion, for square well potential
# Gil-Villegas eqn (34)
def a1SW(n):
    return a1VDW(n)*gHS_eff(n)

# Derivate of a1SW wrt eta
def da1SW_deta(n):

    # derivate of a1VDW wrt eta
    da1VDW_deta = -4*epsilon*(lambdaSW**3-1)

    return gHS_eff(n)*da1VDW_deta + a1VDW(n)*dgHS_eff_deta(n)

# Hard sphere correlation function, effective
# Gil-Villegas eqn (33)
def gHS_eff(n):
    return (1-0.5*eta_eff(n))*(1-eta_eff(n))**-3

# derivative of gHS_eff wrt eta
def dgHS_eff_deta(n):
    return dgHS_eff_deta_eff(n)*deta_eff_deta(n)

# derivative of gHS_eff wrt eta_eff
def dgHS_eff_deta_eff(n):
    return -0.5/(1-eta_eff(n))**3 + 3*(1-0.5*eta_eff(n))/(1-eta_eff(n))**4

# van der Waals attractive parameter
# Gil-Villegas eqn (35)
def a1VDW(n):
    return -4*eta(n)*epsilon*(lambdaSW**3-1)

# 2nd term in High-temp perturbation expansion, for square well potential
# Gil-Villegas eqn (38)
def a2(n):
    return 0.5*epsilon*K(n)*eta(n)*da1SW_deta(n)

# Isothermal compressibility
# Gil-Villegas eqn (22)
def K(n):
    return (1-eta(n))**4/(1+4*eta(n)+4*eta(n)**2)

# Ideal gas
# Gil-Villegas eqn (9)
# NB: In Gil-Villegas, there is a term involving the thermal DeBroglie wavelength (Lambda)
## I ignore this here, because Lambda does not have any contribution from density.
## It only serves to add a constant to the free energy; we absorb this later into the chemical potential
def SWfid(T,n):
    return n*k_B*T*(np.log(n) - 1)

# Hard-sphere free energy
# based on Hughes
def SWfhs(T,n):
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

##########################################

# References

# (1)
## 'Application of a renormalization-group threatment to SAFT-VR'
## E. Forte, F. Llovell, L. F. Vegas, J. P. Trusler, A. Galindo
## Journal of Chemical Physics 134, 154102 (2011)

# (2)
## 'Generalized approach to global renormalization-group theory for fluids'
## A. Sai Venkata Ramana, S. V. G. Menon
## Phys Rev E 85, 041108 (2012)

# (3)
## 'A classical density-functional theory for describing water interfaces'
## J. Hughes, E. J. Krebs, D. Roundy
## J. Chem. Phys. 138, 024509 (2013)

# (4)
## 'SAFT for chain molecules with attractive potentials of variable range'
## A. Gil-Villegas, A. Galindo, P. J. Whitehead, S. J. Mills and G. Jackson
## J. Chem. Phys. 106, 4168 (1997)
