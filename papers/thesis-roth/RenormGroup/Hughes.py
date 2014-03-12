from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Sept 2013

# From "A Classical Density-Functional Theory for Describing Water Interfaces"
# Jessica Hughes, Eric J. Krebs & David Roundy
# J. Chem. Phys. 138, 024509 (2013)

# SAFT Helmholtz free energy is sum of four terms
# Each is defined in Section II of the paper
# This code considers the homogeneous case (density independent of position)
# each term here gives free energy per volume (i.e. f = F/V)

# For water at STP, n = 3.34e22 atoms cm^-3
# or, in atom bohr^-3:
nw = 4.9388942e-3 # (atom bohr^-3

# Constants
# angstrom = 1e-8 # one angstrom in cm
angstrom = 1.8897261 # An angstrom in bohr
Ha = 4.395e-18 # Hartree in Joules
atm = 101325 # atm in Pascal
# k_B = 1.38e-16 # Boltzmann constant, erg K^-1
k_B = 3.16681539628059e-6 # Boltzmann's constant in Hartree/Kelvin
sigma = 3.0342*angstrom # hard-sphere diameter
R = sigma/2 # hard-sphere radius
m = 18.01528 # molar mass of water, g mol^-1
bohr = 5.29177208e-9 # bohr radius in cm
N_A = 6.02213e23 #  Avagadro's number (mol^-1)

# kap_a and eps_a come from "Developing optimal Wertheim-like modesl of water for ue in SAFT and related approaches"
# G. N. I. Clark et al.
# Molecular Physics 104, 3561 (2006)

# mutual volume of interaction
kap_a = 1.06673*angstrom**3 # value has dimensions of length^3

# energy of interaction if spheres hydrogen-bond
eps_a = 1400*k_B # (kelvin x Boltzmann)

# Ideal Gas
def fid(T,n):
    # constants
    hbar = 1.05*10**(-27) # cm^2 g s^-1

    me = 9.109e-28 # electron mass in g
    time_atomic = 2.418e-17 # hbar/E_hartree in s

    hbar_atomic = hbar/bohr**2/me*time_atomic # hbar in atomic units
    m_atomic = m/me/N_A # m in atomic units

    # thermal wavelength
    Lambd = ((2*np.pi*hbar_atomic**2)/(m_atomic*k_B*T))**0.5

    # free energy
    return k_B*T*n*(np.log(n*Lambd**3)-1)

# Fundamental measure densities
def n0(n):
    return n

def n1(n):
    return n*R

def n2(n):
    return n*R**2*4*np.pi

def n3(n):
    return n*R**3*(4/3)*np.pi

# Hard-sphere repulsion
def fhs(T,n):
    # Phi 1
    Phi1 = -n0(n)*np.log(1-n3(n))

    # Phi 2
    Phi2 = (n1(n)*n2(n))/(1-n3(n))

    # Phi 3
    Phi3 = n2(n)**3*((n3(n)+(1-n3(n))**2*np.log(1-n3(n)))/(36*np.pi*n3(n)**2*(1-n3(n))**2))

    # free energy
    return k_B*T*(Phi1 + Phi2 + Phi3)

# Dispersion and Association use a Square Well potential, parameterized by an energy and length scale
# Values from "Developing optimal Wertheim-like modesl of water for ue in SAFT and related approaches"
# G. N. I. Clark et al.
# Molecular Physics 104, 3561 (2006)

eps = 250*k_B # energy
lambd = 1.78890 # dimensionless length

# eta is the packing fraction; volume times number density
# In the homogeneous case, eta is the same as n3(n)
def eta(n):
    return n3(n)

# Coefficients for eta_eff, and derivatives
c1 = 2.25855-1.50349*lambd+0.249434*lambd**2
dc1 = -1.50349+2*0.249434*lambd
c2 = -0.669270+1.40049*lambd-0.827739*lambd**2
dc2 = 1.40049-2*0.827739*lambd
c3 = 10.1576-15.0427*lambd+5.30827*lambd**2
dc3 = -15.0427+2*5.30827*lambd

# effective packing fraction
def eta_eff(n):
    return c1*eta(n) + c2*eta(n)**2 + c3*eta(n)**3

# derivative of eta_eff wrt eta
def deta_eff_deta(n):
    return c1+2*c2*eta(n)+3*c3*eta(n)**2

# derivative of eta_eff wrt lambda
def deta_eff_dlambd(n):
    return dc1*eta(n) + dc2*eta(n)**2 + dc3*eta(n)**3

# van der Waals attractive parameter
def a1VDW(n):
    return -4*eta(n)*eps*(lambd**3-1)

# derivate of a1VDW wrt eta
da1VDW_deta = -4*eps*(lambd**3-1)

# derivative of a1VDW wrt lambda
def da1VDW_dlambd(n):
    return -12*eps*lambd**2*eta(n)

# Hard sphere correlation function, effective
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

# Isothermal compressibility
def K(n):
    return (1-eta(n))**4/(1+4*eta(n)+4*eta(n)**2)

# 1st term in High-temp perturbation expansion, for square well potential
def a1SW(n):
    return a1VDW(n)*gHS_eff(n)

# Derivate of a1SW wrt eta
def da1SW_deta(n):
    return gHS_eff(n)*da1VDW_deta + a1VDW(n)*dgHS_eff_deta(n)

# Derivative of a1SW wrt lambda
def da1SW_dlambd(n):
    return da1VDW_dlambd(n)*gHS_eff(n) + a1VDW(n)*dgHS_eff_dlambd(n)

# 2nd term
def a2(n):
    return 0.5*eps*K(n)*eta(n)*da1SW_deta(n)

# Dispesion Free Energy
def fdisp(T,n):
    return (a1SW(n) + 1/(k_B*T)*a2(n))*n

# correlation functions
def gHS(n):
    return 1/(1-n3(n))+ (R/2)*n2(n)/(1-n3(n))**2 + (R**2/18)*n2(n)**2/(1-n3(n))**3

def gSW(T,n):
    return gHS(n) + 1/(4*k_B*T)*(da1SW_deta(n) - lambd/(3*eta(n))*da1SW_dlambd(n))

# Delta is a measure of hydrogen-bonding proability
def Delta(T,n):
    return kap_a*gSW(T,n)*(np.exp(eps_a/(k_B*T)) - 1)

# Fraction of association sites NOT hydrogen-bonded
def X(T,n):
    return ((1 + 8*n0(n)*Delta(T,n))**0.5 - 1)/(4*n0(n)*Delta(T,n))

#  Association Free Energy
def fassoc(T,n):
    return 4*k_B*T*n0(n)*(np.log(X(T,n)) - X(T,n)/2 + 0.5)

# Total free energy
def f(T,n):
    return fid(T,n)+fhs(T,n)+fdisp(T,n)+fassoc(T,n)

# Derivative wrt number at const temp
def df_dn(T,n):
    # step size
    dn = 1e-8*n
    return (f(T,n + dn/2) - f(T,n - dn/2))/dn

# Pressure
def p(T,n):
    return n*df_dn(T,n) - f(T,n)

# Conversion factor taking n from atoms/bohr^3 to g/mL
conv_n = m/bohr**3/N_A # g cm^-3 atoms^-1

# Conversion factor taking p from Hartree/bohr^3 to atm
conv_p = Ha*1e6/bohr**3/atm

# Grand free energy
def Phi(T,n,prefac):
    '''Grand free energy'''
    nl = prefac*nw # we will adjust the prefactor to find accurate chemical potential
    mu = df_dn(T,nl) # chemical potential; derivative of free energy around a particular density
    return f(T,n) - mu*n
