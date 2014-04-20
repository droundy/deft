from __future__ import division
import RG
import sys
import minmax_RG
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

T = 0.1
nparticular = 0.025
N = 20 # For publication plots, make this bigger (40? 80? 100? You decide!)
Tc = 1.33
Tlow = T

# Recursion depth
iterations = 0

# Bounds of minimization.
# Use range that works for many temperatures
# Left (vapor)
a_vap = 1e-10/(RG.sigma**3*np.pi/6)
c_vap = nparticular
# Right (liquid)
a_liq = nparticular
c_liq = 0.55/(RG.sigma**3*np.pi/6)
###############################

# Open file for output
fout = open('figs/npart_RG-out.dat','w')

# label the columns of the output
fout.write('#T     nvapor     nliquid       phi(nvap)        phi(nliq)         nparticular\n')#       phi_avg')

# Do first temperature before the loop
nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,nparticular,iterations)
print 'nvap,phi_vap',nvapor,phi_vapor
nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,nparticular,iterations)
print 'nliq,phi_liq',nliquid,phi_liquid
sys.stdout.flush()


#while T < 1.4:
for j in xrange(0,N+1):
    print ' Top of for loop; j =',j
    T = (Tc - Tlow)*(1 - ((N-j)/N)**4) + Tlow

    fout.flush()
    # Starting point for new nparticular is abscissa of max RG.phi with old nparticular
    nparticular = minmax_RG.maximize(RG.phi,T,nvapor, nliquid, nparticular,iterations)
    print '   npart =',nparticular

    # I'm looking at the minima of RG.phi
    c_vap = nparticular
    a_liq = nparticular

    nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,c_vap,a_vap,nparticular,iterations)
    nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,nparticular,iterations)

    tol = 1e-5
    fnpart = RG.phi(T, nparticular, nparticular, iterations)

    # Compare the two minima in RG.phi
    while np.fabs(phi_vapor - phi_liquid)/np.fabs(fnpart) > tol:
        print '     in the loop comparing two minima in RG.phi'

        delta_mu = (phi_liquid - phi_vapor)/(nliquid - nvapor)
        if delta_mu > 0:
            print '     delta_mu is positive'
        elif delta_mu < 0:
            print '     delta_mu is negative'
        else:
            print '     delta_mu is zero'

        def newphi(T, n, npart, i):
            return RG.phi(T, n, npart, i) - delta_mu*n

        print '     finding max in newphi'
        nparticular = minmax_RG.maximize(newphi,T,nvapor,nliquid,nparticular,iterations)
        print '     new nparticular is',nparticular

        fnpart = RG.phi(T, nparticular, nparticular, iterations)
        print '     fnpart =',fnpart

        print '     finding new nvapor,phi_vapor'
        nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,nparticular,iterations)

        print '     finding new nliquid,phi_liquid'
        nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,nparticular,iterations)

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
    print '   T, etaVap, etaLiq',T,nvapor/(RG.sigma**3*np.pi/6),nliquid/(RG.sigma**3*np.pi/6)
    print ' Going back to top of for loop'
