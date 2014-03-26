from __future__ import division
import SW
import sys
import minmax
import numpy as np
import pylab

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Dec 2013

# General scheme:
# 1) Find nparticular at an easy temp (say, 300k) by hand
# ___
# |2) Using that nparticular, find maximum in Phi
# |3) The density at the maximum (which should be similar to the nparticular used)
# |   is your first guess for the nparticular at a slightly higher temperature
# |4) Starting from that guess, adjust the nparticular until the two minima in Phi are equal
# |    4a) npart inc. --> right-side minimum inc.
# |    4b) npart dec. --> right-side minimum dec.
# |    The left-side minimum seems to not change as drastically
# |5) go back to step (2)
# |__

###############################
# Initial conditions; dependent on you system
# T = 0.001 # SW units
# nparticular = 0.14566 # I found this by hand

T = 0.1
nparticular = 0.05/(SW.sigma**3*np.pi/6)
N = 20 # For publication plots, make this bigger (40? 80? 100? You decide!)
Tc = 1.33
Tlow = T

# Bounds of minimization.
# Use range that works for many temperatures
# Left (vapor)
a_vap = 1e-10/(SW.sigma**3*np.pi/6)
c_vap = nparticular
# Right (liquid)
a_liq = nparticular
c_liq = 0.55/(SW.sigma**3*np.pi/6)
###############################

# Open file for output
fout = open('figs/npart_SW-out.dat','w')

# label the columns of the output
fout.write('#T     nvapor     nliquid       phi(nvap)        phi(nliq)         nparticular\n')#       phi_avg')

# Do first temperature before the loop
nvapor,phi_vapor = minmax.minimize(SW.phi,T,a_vap,c_vap,nparticular)
nliquid,phi_liquid = minmax.minimize(SW.phi,T,a_liq,c_liq,nparticular)
sys.stdout.flush()

#while T < 1.4:
for i in xrange(0,N+1):
    T = (Tc - Tlow)*(1 - ((N-i)/N)**4) + Tlow

    fout.flush()
    # Starting point for new nparticular is abscissa of max SW.phi with old nparticular
    nparticular = minmax.maximize(SW.phi,T,nvapor, nliquid, nparticular)

    # I'm looking at the minima of SW.phi
    c_vap = nparticular
    a_liq = nparticular

    nvapor,phi_vapor = minmax.minimize(SW.phi,T,c_vap,a_vap,nparticular)
    nliquid,phi_liquid = minmax.minimize(SW.phi,T,a_liq,c_liq,nparticular)

    tol = 1e-5
    fnpart = SW.phi(T, nparticular, nparticular)

    # Compare the two minima in SW.phi
    while np.fabs(phi_vapor - phi_liquid)/np.fabs(fnpart) > tol:
        delta_mu = (phi_liquid - phi_vapor)/(nliquid - nvapor)
        def newphi(T, n, npart):
            return SW.phi(T, n, npart) - delta_mu*n
        nparticular = minmax.maximize(newphi,T,nvapor,nliquid,nparticular)
        fnpart = SW.phi(T, nparticular, nparticular)

        nvapor,phi_vapor = minmax.minimize(SW.phi,T,a_vap,c_vap,nparticular)
        nliquid,phi_liquid = minmax.minimize(SW.phi,T,a_liq,c_liq,nparticular)

    fout.write(str(T))
    fout.write('  ')
    fout.write(str(nvapor))
    fout.write('  ')
    fout.write(str(nliquid))
    fout.write('  ')
    fout.write(str(phi_vapor))
    fout.write('  ')
    fout.write(str(phi_liquid))
    fout.write('  ')
    fout.write(str(nparticular))
    fout.write('\n')
    sys.stdout.flush();
    print T,'\t',nvapor/(SW.sigma**3*np.pi/6),'\t',nliquid/(SW.sigma**3*np.pi/6)

    # Set temp slightly higher
#    T += 0.01
