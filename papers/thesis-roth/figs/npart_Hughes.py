from __future__ import division
import Hughes as H
import minmax
import numpy as np

"""Find the optimal 'prefac' to get the chemical potenential for coexistance in the function Hughes.Phi"""

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Dec 2013

# General scheme:
# 1) Find prefac at an easy temp (say, 300k) by hand
# ___
# |2) Using that prefac, find maximum in H.Phi
# |3) The density at the maximum (which should be similar to the prefac used)
# |   is your first guess for the prefac at a slightly higher temperature
# |4) Starting from that guess, adjust the prefac until the two minima in H.Phi are equal
# |    4a) prefac inc. --> right-side minimum inc.
# |    4b) prefac dec. --> right-side minimum dec.
# |    The left-side minimum seems to not change as drastically
# |5) go back to step (2)
# |__

###############################
# Initial conditions:
T = 276 # Kelvin
nparticular = 0.0022182 # I found this by hand; initial prefac from above times n in atom/bohr^3 at STP


N = 20 # For publication plots, make this bigger (40? 80? 100? You decide!)
Tc = 700
Tlow = T

# IC's for bounds of minimization.
# Use range that works for many temperatures
# Left
a_vap = 1e-8/H.conv_n
c_vap = 0.4/H.conv_n
# Right
a_liq = 0.4/H.conv_n
c_liq = 1.05/H.conv_n
###############################


# Open file for output
fout = open('figs/npart_Hughes-out.dat','w')

# label the columns of the output
fout.write('#T    nvapor(g/ml)    nliquid(g/mL)    phi(nvap)(atm)    phi(nliq)(atm)    nparticular')

# Do first temperature before the loop
leftmin_n,leftmin_phi = minmax.minimize(H.Phi_alt,T,a_vap,c_vap,nparticular)
rightmin_n,rightmin_phi = minmax.minimize(H.Phi_alt,T,a_liq,c_liq,nparticular)


leftmin_n = a_vap
rightmin_n = c_liq

# while T < 800:
for i in xrange(0,N+1):
    T = (Tc - Tlow)*(1 - ((N-i)/N)**4) + Tlow

    fout.flush()
    # Starting point for new nparticular is abscissa of max H.Phi_alt with old nparticular
    nparticular = minmax.maximize(H.Phi_alt,T,leftmin_n, rightmin_n, nparticular)
    tol_nparticular = 1e-2/H.conv_n


    # I'm looking at the minima of H.Phi_alt
    c_vap = nparticular
    a_liq = nparticular
    leftmin_n,leftmin_phi = minmax.minimize(H.Phi_alt,T,c_vap,a_vap,nparticular)
    rightmin_n,rightmin_phi = minmax.minimize(H.Phi_alt,T,a_liq,c_liq,nparticular)
    tol_min = 0.5/H.conv_p

    # Compare the two minima in H.Phi_alt
    while np.fabs(leftmin_phi - rightmin_phi) > tol_min:
        delta_mu = (rightmin_phi - leftmin_phi)/(rightmin_n - leftmin_n)
        def newphi(T, n, npart):
            return H.Phi_alt(T, n, npart) - delta_mu*n
        oldnparticular = nparticular
        nparticular = minmax.maximize(newphi,T,leftmin_n,rightmin_n,nparticular)

        leftmin_n,leftmin_phi = minmax.minimize(H.Phi_alt,T,a_vap,c_vap,nparticular)
        rightmin_n,rightmin_phi = minmax.minimize(H.Phi_alt,T,a_liq,c_liq,nparticular)

    fout.write('\n')
    fout.write(str(T))
    fout.write('  ')
    fout.write(str(leftmin_n*H.conv_n))
    fout.write('  ')
    fout.write(str(rightmin_n*H.conv_n))
    fout.write('  ')
    fout.write(str(leftmin_phi*H.conv_p))
    fout.write('  ')
    fout.write(str(rightmin_phi*H.conv_p))
    fout.write('  ')
    fout.write(str(nparticular))
    print T

    # Set temp slightly higher
#    dT = 5000/T
#    T = T + dT
