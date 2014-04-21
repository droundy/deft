from __future__ import division
import numpy as np

# Author: Dan Roth
# Email: Daniel.Edward.Roth@gmail.com
# Date: February 2014

# Free energy density for Square Well liquid
# Based on theory presented it:
## SAFT for chain molecules with attractive potentials of variable range
## A. Gil-Villegas, A. Galindo, P. J. Whitehead, S. J. Mills and G. Jackson
## J. Chem. Phys. 106, 4168 (1997)

# I also reference the following:
## A Classical Density-Functional Theory for Describing Water Interfaces
## J. Hughes, E. J. Krebs and D. Roundy
## J. Chem. Phys. 138, 024509 (2013)

# Constants
k_B = 1 # define Boltzman constant as unity, since we will be working in SW units

# Square Well
sigma = 2 # HS diameter
epsilon = 1 # depth of well
lambdaSW = 1.5 # range of well
R = sigma/2 # HS radius

# Coefficients for eta_eff, and derivatives
# Gil-Villegas eqn (27)
c1 = 2.25855-1.50349*lambdaSW+0.249434*lambdaSW**2
dc1 = -1.50349+2*0.249434*lambdaSW
c2 = -0.669270+1.40049*lambdaSW-0.827739*lambdaSW**2
dc2 = 1.40049-2*0.827739*lambdaSW
c3 = 10.1576-15.0427*lambdaSW+5.30827*lambdaSW**2
dc3 = -15.0427+2*5.30827*lambdaSW

# derivate of a1VDW wrt eta
# This is required for fdisp
da1VDW_deta = -4*epsilon*(lambdaSW**3-1)

# Total free energy per volume
# NB: Gil-Villegas also includes a Chain term; we do not deal with chains, so I leave it off
def f(T,n):
    return fid(T,n) + fdisp(T,n) + fhs(T,n)

# Grand free energy per volume
def phi(T,n,nparticular):
#    mu = df_dn(T,nparticular)
    mu = alt_df_dn(T,nparticular)
    return f(T,n) - mu*n

# Derivative of free energy wrt number at const temp
def alt_df_dn(T,n):
    # step size
    dn = 1e-6*n
    return (f(T,n + dn/2) - f(T,n - dn/2))/dn

# Derivative of free energy wrt number (const. temp)
def df_dn(T,n):

    # Constants for eta_eff
    # Redefine here, rather than using SW.c1 (etc.) in case the value of lambdaSW is different between SW.py and this program
    c1 = 2.25855-1.50349*lambdaSW+0.249434*lambdaSW**2
    c2 = -0.669270+1.40049*lambdaSW-0.827739*lambdaSW**2
    c3 = 10.1576-15.0427*lambdaSW+5.30827*lambdaSW**2

    # Sim. for n1,n2,n3
    n0 = n
    n1 = R*n
    n2 = R**2*4*np.pi*n
    n3 = R**3*4/3*np.pi*n

    # n3 compliment
    n3c = 1 - n3

    # These pop up a lot
    dn3_dn = 4/3*np.pi*R**3
    deta_dn = dn3_dn # for homogeneous case deta_dn is the same as dn3_dn
    deta_eff_dn = deta_dn*(c1 + 2*c2*eta(n) + 3*c3*eta(n)**2)

    # ideal gas term
    dfid_dn = k_B*T*np.log(n)

    # Hard sphere term
    # This gets a bit complicated; bear with me
    dPhi1_dn = n0/n3c*dn3_dn - np.log(n3c)
    dPhi2_dn = 1/n3c*(n2*R + n1*4*np.pi*R**2) - n1*n2/n3c**2*dn3_dn
    dPhi3_dn = 1/36/np.pi/n3**2/n3c**2*((n3 + n3c**2*np.log(n3c))*3*n2**2*4*np.pi*R**2 + n2**3*dn3_dn*(1 - 2*n3c*np.log(n3c) - n3*n3c) - (n2**3*(n3 + n3c**2*np.log(n3c)))/n3/n3c*dn3_dn*(2*n3c - n*n3))

    dfhs_dn = k_B*T*(dPhi1_dn + dPhi2_dn + dPhi3_dn)

    # This is a bit more complicated for dfdisp_dn...just hang in there, you'll make it through!

    # Term from derivative a1SW
    da1SW_dn = deta_dn*(a1VDW(n)/(1 - eta_eff(n))**3*(c1 + 2*c2*eta(n) + 3*c3*eta(n)**2)*(3*(1 - 0.5*eta_eff(n))/(1 - eta_eff(n)) - 0.5) - 4*epsilon*(lambdaSW**3 - 1)*gHS_eff(n))
    dK_dn = -(1 - eta(n))**3/(1 + 4*eta(n) + 4*eta(n)**2)*deta_dn*((1 - eta(n))/(1 + 4*eta(n) + 4*eta(n)**2)*(8*eta(n) + 4) + 4)
    da2_dn = -2*epsilon**2*(lambdaSW**3 - 1)*(eta(n)*dK_dn + K(n)*deta_dn)

    dfdisp_dn = a1SW(n) + 1/k_B/T*a2(n) + n*(da1SW_dn + 1/k_B/T*da2_dn)

    return dfid_dn + dfdisp_dn + dfhs_dn

# Ideal gas
# Gil-Villegas eqn (9)
# NB: In Gil-Villegas, there is a term involving the thermal DeBroglie wavelength (Lambda)
## I ignore this here, because Lambda does not have any contribution from density.
## It only serves to add a constant to the free energy; we absorb this later into the chemical potential
def fid(T,n):
    return n*k_B*T*(np.log(n) - 1)

# Hard-sphere free energy
# based on Hughes
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

# Dispersion
# Hughes calls this Dispersion, Gil-Villegas calls it Monomer
# This is taken from the code based on Hughes (Hughes.py)
def fdisp(T,n):
    return (a1SW(n) + 1/(k_B*T)*a2(n))*n

# 1st term in High-temp perturbation expansion, for square well potential
# Gil-Villegas eqn (34)
def a1SW(n):
    return a1VDW(n)*gHS_eff(n)

# Derivate of a1SW wrt eta
def da1SW_deta(n):
    return gHS_eff(n)*da1VDW_deta + a1VDW(n)*dgHS_eff_deta(n)

# Derivative of a1SW wrt lambda
def da1SW_dlambd(n):
    return da1VDW_dlambd(n)*gHS_eff(n) + a1VDW(n)*dgHS_eff_dlambd(n)

# derivative of a1VDW wrt lambda
def da1VDW_dlambd(n):
    return -12*epsilon*lambdaSW**2*eta(n)

# 2nd term in High-temp perturbation expansion, for square well potential
# Gil-Villegas eqn (38)
def a2(n):
    return 0.5*epsilon*K(n)*eta(n)*da1SW_deta(n)

# Hard sphere correlation function, effective
# Gil-Villegas eqn (33)
def gHS_eff(n):
    return (1-0.5*eta_eff(n))*(1-eta_eff(n))**-3

# derivative of gHS_eff wrt eta_eff
def dgHS_eff_deta_eff(n):
    return -0.5/(1-eta_eff(n))**3 + 3*(1-0.5*eta_eff(n))/(1-eta_eff(n))**4

# derivative of gHS_eff wrt eta
def dgHS_eff_deta(n):
    return dgHS_eff_deta_eff(n)*deta_eff_deta(n)

# derivative of gHS_eff wrt lamda
def dgHS_eff_dlambd(n):
    return dgHS_eff_deta_eff(n)*deta_eff_dlambd(n)

# effective packing fraction
# Gil-Villegas eqn (36)
def eta_eff(n):
    return c1*eta(n) + c2*eta(n)**2 + c3*eta(n)**3

# derivative of eta_eff wrt eta
def deta_eff_deta(n):
    return c1+2*c2*eta(n)+3*c3*eta(n)**2

# derivative of eta_eff wrt lambda
def deta_eff_dlambd(n):
    return dc1*eta(n) + dc2*eta(n)**2 + dc3*eta(n)**3

# eta is the packing fraction; volume times number density
# In the homogeneous case, eta is the same as n3(n) (from FMT; see Hughes 2013)
def eta(n):
    return n*R**3*(4/3)*np.pi

# van der Waals attractive parameter
# Gil-Villegas eqn (35)
def a1VDW(n):
    return -4*eta(n)*epsilon*(lambdaSW**3-1)

# Isothermal compressibility
# Gil-Villegas eqn (22)
def K(n):
    return (1-eta(n))**4/(1+4*eta(n)+4*eta(n)**2)

# Pressure, defined in Hughes' thesis, eqn 38
def P(T,n):
    return n*df_dn(T,n) - f(T,n)
