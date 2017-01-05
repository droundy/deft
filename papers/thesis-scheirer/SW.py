from __future__ import division
import numpy as np

###############################################################################################
# Author: Ryan Scheirer                                                                       #
# Email: scheirer@oregonstate.edu                                                             #
# Date: February 2016                                                                         #
                                                                                              #
# Uses SAFT to calculate Free energy density for Square Well liquid based off of the paper... #
# ... SAFT for chain molecules with attractive potentials of variable range...                #
# ... A. Gil-Villegas, A. Galindo, P.J. Whitehead, S. J. Mills and G. Jackson...              #
# ... J. Chem. Phys. 106, 4168 (1997)                                                         # 
                                                                                              #
# This program also references the following paper...                                         # 
# ... A Classical Density-Functional Theory for Describing Water Interfaces...                #
# ... J. Hughes, E. J. Krebs and D. Roundy...                                                 #
# ... J. Chem. Phys. 138, 024509 (2013)                                                       #
#                                                                                             #
# This program largley uses Dan Roth's SW.py code, I simply edit this accordingly.            #    
#                                                                                             #
# Also included is the necessary constants and definitions to calculate the...                #
# ...Van Der Waal's free energy density.                                                      # 
###############################################################################################





################################## START CONSTANTS ############################################
#                                                                                             #
#                                                                                             #
kb = 1          # define Boltzman constant as unity, since we will be working in SW units


### Square Well Stuff ###
sigma = 2       # HS diameter
epsilon = 1     # depth of well
lambdaSW = 1.5  # range of interaction (valid range: 1.3 - 1.8) (NOTE: sigma*lambdaSW - sigma = well width)
R = sigma/2     # HS radius

# Coefficients for eta_eff (effective filling fraction) A. Gil-Villegas eqn(37)
c1 = 2.25855 - 1.50349*lambdaSW + 0.249434*lambdaSW**2
c2 = -0.669270 + 1.40049*lambdaSW - 0.827739*lambdaSW**2
c3 = 10.1576 - 15.0427*lambdaSW + 5.30827*lambdaSW**2
# Coefficients for d(eta_eff)_d(lambaSW) (derivatives of the above 3 coefficients w.r.t. lambdaSW)
dc1 = -1.50349 + 2*0.249434*lambdaSW
dc2 = 1.40049 - 2*0.827739*lambdaSW
dc3 = -15.0427 + 2*5.30827*lambdaSW

# derivate of a1VDW wrt eta (this is required for fdisp) A. Gil-Villegas eqn(35)
da1VDW_deta = -4*epsilon*(lambdaSW**3-1)


### Van Der Waal's Stuff ###
# Only used these to test SnAFT against Van Der Waal's but due to time limitations I never fully explored this fun comparison
B=(np.pi*sigma**3)/6.0
alphaVDW = ((2*np.pi*epsilon*sigma**3)/3)*(lambdaSW**3-1)
#                                                                                             #
#                                                                                             #
################################### END CONSTANTS #############################################





################################### START FREE ENERGIES #######################################
#                                                                                             #
#                                                                                             #
def ftot(T,n):
    ### Totl free energy per volume ### A. Gil-Villegas eqn(3)
    # A. Gil-Villegas also includes Chain and an Association terms; we do not deal with these, thus I leave them out
    # A. Gil-Villegas also uses mono where we break that up into HS and disp (so mono = HS + disp)
    # n is the number density
    # T is reduced temperature
    return fid(T,n) + fhs(T,n) + fdisp(T,n)


def frg(T,n):
    # Was using this as a test function when comparing my RG code.
    # This is considered f0, SnAFT minus one of the perturbation terms found in the dispersion free energy. 
    return fid(T,n) + fhs(T,n) + 1/(kb*T)*a2(n)*n


def fxs(T,n):
    ### Excess free enery per volume (used for initial comparisons with M.C. simulations) ###
    return fdisp(T,n)


def fxsc(T,n):
    ### Excess free energy per volume corrected so that at T=infinity the slope is zero (used for initial comparisons with M.C. simulations) ###
    return fxs(T,n) - num_dfxs_dT(100,n)*T


def num_dftot_dT(T,n):
    ### Derivative of total free energy per volume w.r.t. T  (NOTE: at fixed volume and pressure, this is neg. entropy) ###
    dT=(0.001)*np.abs(n)
    df1=ftot(T+dT,n)-ftot(T-dT,n-dn)
    df2=ftot(T+2*dT,n)-ftot(T-2*dT,n)
    df3=ftot(T+3*dT,n)-ftot(T-3*dT,n)
    df4=ftot(T+4*dT,n)-ftot(T-4*dT,n)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dT


def num_dfxs_dT(T,n):
    ### Derivative of excess free energy per volume w.r.t. T  (NOTE: at fixed volume and pressure, this is entropy_ex) ###
    dT=(0.001)*np.abs(n)
    df1=fxs(T+dT,n)-fxs(T-dT,n-dn)
    df2=fxs(T+2*dT,n)-fxs(T-2*dT,n)
    df3=fxs(T+3*dT,n)-fxs(T-3*dT,n)
    df4=fxs(T+4*dT,n)-fxs(T-4*dT,n)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dT


def numL_dftot_dn(T,n):
    ### Derivative of total free energy per volume w.r.t. n  (NOTE: at fixed volume and pressure, this is mu) ###
    dn = 1e-5*n
    return (ftot(T,n + dn) - ftot(T,n - dn))/(2*dn)
    

def numH_dftot_dn(T,n):
    ### Derivative of total free energy per volume w.r.t. n  (NOTE: at fixed volume and pressure, this is mu) ###
    # This higher order derivative was tested against, and works better than numL_dftot_dn
    dn=(0.001)*np.abs(n)
    df1=ftot(T,n+dn)-ftot(T,n-dn)
    df2=ftot(T,n+2*dn)-ftot(T,n-2*dn)
    df3=ftot(T,n+3*dn)-ftot(T,n-3*dn)
    df4=ftot(T,n+4*dn)-ftot(T,n-4*dn)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dn


def num_dfxs_dn(T,n):
    ### Derivative of excess free energy per volume w.r.t. n  (NOTE: at fixed volume and pressure, this is mu_ex) ###
    dn=(0.001)*np.abs(n)
    df1=fxs(T,n+dn)-fxs(T,n-dn)
    df2=fxs(T,n+2*dn)-fxs(T,n-2*dn)
    df3=fxst(T,n+3*dn)-fxs(T,n-3*dn)
    df4=fxs(T,n+4*dn)-fxs(T,n-4*dn)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dn


def phi(T,n,npart):
    ### Grand free energy per volume ###
    mu = numH_dftot_dn(T,npart)
    return ftot(T,n)-mu*n
    

def fvdw(T,n):
    ### Van Der Waal's free energy per volume ###
    return fid(T,n) - n*kb*T*(np.log(1-4*n*B))-(n**2)*alphaVDW

	
def numH_dfvdw_dn(T,n):
    ### Derivative of Van Der Waal's free energy per volume w.r.t. n ###
    dn=(0.001)*np.abs(n)
    df1=fvdw(T,n+dn)-fvdw(T,n-dn)
    df2=fvdw(T,n+2*dn)-fvdw(T,n-2*dn)
    df3=fvdw(T,n+3*dn)-fvdw(T,n-3*dn)
    df4=fvdw(T,n+4*dn)-fvdw(T,n-4*dn)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dn


def fmvdw(T,n):
    ### Modified Van Der Waal's free energy per volume ###
    return fid(T,n) + n*kb*T*((4*n*B-3*(n**2)*(B**2))/((1-n*B)**2)) - (n**2)*alphaVDW

	
def numH_dfmvdw_dn(T,n):
    ### Derivative of modified Van Der Waal's free energy per volume w.r.t. n ###
    dn=(0.001)*np.abs(n)
    df1=fmvdw(T,n+dn)-fmvdw(T,n-dn)
    df2=fmvdw(T,n+2*dn)-fmvdw(T,n-2*dn)
    df3=fmvdw(T,n+3*dn)-fmvdw(T,n-3*dn)
    df4=fmvdw(T,n+4*dn)-fmvdw(T,n-4*dn)
    return (df1*4.0/5.0-df2/5.0+df3*4.0/105.0-df4/280.0)/dn

#                                                                                             #
#                                                                                             #
###################################### END FREE ENERGIES ######################################




    
################################## START THERMODYNAMIC PROPERTIES #############################
#                                                                                             #
#                                                                                             #
def findSxs(T,n):
    ### Excess free energy (used to compare with M.C. simulations)
    return -num_dfxs_dT(T,n)

def findStot(T,n):
    ### Total entropy ###
    return -num_dftot_dT(T,n)
    
    
def findP(T,n):
    ### Pressure ###
    return n*numH_dftot_dn(T,n) - ftot(T,n)
#                                                                                             #
#                                                                                             #
################################# END THERMODYNAMIC PROPERTIES ################################





###################################### START IDEAL GAS ########################################
#                                                                                             #
#                                                                                             #
def fid(T,n):
    # A. Gil-Villegas eqn(9)
    # Naturally we set lambda to 1 of course!
    return n*kb*T*(np.log(n) - 1)
#                                                                                             #
#                                                                                             #
######################################### END IDEAL GAS #######################################





######################################## START HARD-SPHERE ####################################
#                               NOTE: definitions are based on J. Hughes et al.               #
#                                                                                             #
def fhs(T,n):
    # Phi 1: J. Hughes eqn(7)
    Phi1 = -n0(n)*np.log(1-n3(n))

    # Phi 2: J. Hughes eqn(8)
    Phi2 = (n1(n)*n2(n))/(1-n3(n))

    # Phi 3: J. Hughes eqn(9)
    Phi3 = n2(n)**3*((n3(n)+(1-n3(n))**2*np.log(1-n3(n)))/(36*np.pi*n3(n)**2*(1-n3(n))**2))

    # HS free energy
    return kb*T*(Phi1 + Phi2 + Phi3)


# Fundamental measure densities
def n0(n):
    return n

def n1(n):
    return n*R

def n2(n):
    return n*R**2*4*np.pi

def n3(n):
    return n*R**3*(4/3)*np.pi
#                                                                                             #
#                                                                                             #
###################################### END HARD-SPHERE ########################################





###################################### START DISPERSION #######################################
#                                                                                             #
#                                                                                             #
def fdisp(T,n):
    # A. Gil-Villegas, basically eqn(18) and eqn(10) combined
    return (a1SW(n) + 1/(kb*T)*a2(n))*n



### a1SW(n) and a2(n) common stuff
def a1VDW(n):
    # van der Waals attractive parameter
    # A. Gil-Villegas eqn(35)
    return -4*eta(n)*epsilon*(lambdaSW**3-1)

def gHS_eff(n):
    # Hard sphere correlation function, effective
    # A. Gil-Villegas eqn(33)
    return (1-0.5*eta_eff(n))*(1-eta_eff(n))**-3

def eta(n):
    # eta is the packing fraction; volume times number density
    # In the homogeneous case, eta is the same as n3(n) (from FMT; see Hughes 2013)
    return n*R**3*(4/3)*np.pi

def eta_eff(n):
    # effective packing fraction
    # A. Gil-Villegas eqn(36)
    return c1*eta(n) + c2*eta(n)**2 + c3*eta(n)**3




### a1SW(n) stuff
def a1SW(n):
    # 1st term in High-temp perturbation expansion, for square well potential
    # A. Gil-Villegas eqn(34)
    return a1VDW(n)*gHS_eff(n)




def deta_eff_deta(n):
    # derivative of eta_eff wrt eta
    return c1+2*c2*eta(n)+3*c3*eta(n)**2





### a2(n) stuff
def a2(n):
    # 2nd term in High-temp perturbation expansion, for square well potential
    # A. Gil-Villegas eqn(38)
    return 0.5*epsilon*K(n)*eta(n)*da1SW_deta(n)

def K(n):
    # Isothermal compressibility
    # A. Gil-Villegas eqn(22)
    return (1-eta(n))**4/(1+4*eta(n)+4*eta(n)**2)

def da1SW_deta(n):
    # Derivate of a1SW w.r.t. eta
    return gHS_eff(n)*da1VDW_deta + a1VDW(n)*dgHS_eff_deta(n)

def dgHS_eff_deta(n):
    # derivative of gHS_eff wrt eta
    return dgHS_eff_deta_eff(n)*deta_eff_deta(n)

def dgHS_eff_deta_eff(n):
    # derivative of gHS_eff w.r.t. eta_eff
    return -0.5/(1-eta_eff(n))**3 + 3*(1-0.5*eta_eff(n))/(1-eta_eff(n))**4



#                                                                                             #
#                                                                                             #
####################################### END DISPERSION ########################################
