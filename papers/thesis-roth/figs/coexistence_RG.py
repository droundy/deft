from __future__ import division
import RG_v2 as RG
import minmax_RG_v2 as minmax_RG
import numpy as np

# This needs to be so for SCons
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')

import matplotlib.pyplot as plt

# The following is not from the standard library. Find it in:
# (...)/deft/papers/thesis-roth/figs
import forte

##############################################################################################
# Author: Dan Roth                                                                           #
# E-mail: Daniel.Edward.Roth@gmail.com                                                       #
# Date: June 2014                                                                            #
#                                                                                            #
# This is the second-ish version (previously called npart_RG)                                #
# I have changed the main RG program, so this needs to be updated accordingly                #
#                                                                                            #
# General scheme:                                                                            #
# 1) Find a value for mu by hand                                                             #
# ___                                                                                        #
# |2) The previous mu value is your first guess for mu at a slightly higher temperature      #
# |3) Starting from that guess, adjust the value of mu:                                      #
# |    3a) mu more negative --> phi(nliquid) increases                                       #
# |    3b) mu less negative --> phi(nliquid) decreases                                       #
# |    The phi(nvapor) does not seem to change as drastically                                #
# |4) go back to step (2)                                                                    #
# |__                                                                                        #
##############################################################################################

# input iteration depth as a parameter
def npart(iterations):
    # Initial conditions; dependent on you system

    T = 0.5
    mu = -5.751 # I found this by plotting. It is good for T=0.5; i=0

    N = 20 # For publication plots, make this bigger (40? 80? 100? You decide!)
    Tc = 1.33 # Rough estimate based on previous plots
    Tlow = T

    # Bounds of minimization.
    # Use range that works for many temperatures
    # Outside bounds
    a_vap = 1e-10/(RG.sigma**3*np.pi/6)
    c_liq = 0.75/(RG.sigma**3*np.pi/6)
    # Use the central max as the inside bound
    nmid = minmax_RG.maximize(RG.phi,T,a_vap,c_liq,mu,iterations)
    c_vap = nmid
    a_liq = nmid

    ###############################

    # Open file for output
    fout = open('figs/coexistence_RG-i%d-out.dat'%iterations,'w')

    # label the columns of the output
    fout.write('#T     nvapor    nmid     nliquid       phi(nvap)       phi(nmid)         phi(nliq)         mu\n')
    sys.stdout.flush()

    # Do first temperature before the loop
    nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,mu,iterations)
    print '  initial nvap,phi_vap',nvapor,phi_vapor
    nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,mu,iterations)
    print '  initial nliq,phi_liq',nliquid,phi_liquid
    print '  initial mu',mu

    #while T < 1.4:
    for j in xrange(0,N+1):
        # Temp points get closer as we near the critical point
        T = (Tc - Tlow)*(1 - ((N-j)/N)**4) + Tlow

        fout.flush()
        # Starting point for new nparticular is abscissa of max RG.phi with old nparticular
        nmid = minmax_RG.maximize(RG.phi,T,nvapor, nliquid, mu,iterations)

        # I'm looking at the minima of RG.phi
        c_vap = nmid
        a_liq = nmid

        nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,c_vap,a_vap,mu,iterations)
        nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,mu,iterations)

        tol = 1e-5
        phi_mid = RG.phi(T, nmid, mu, iterations)

        # Compare the two minima in RG.phi
        print '    entering while loop'
        while np.fabs(phi_vapor - phi_liquid)/np.fabs(phi_mid) > tol:

            delta_mu = (phi_liquid - phi_vapor)/(nliquid - nvapor)
            print '      delta_mu=',delta_mu

            def newphi(T, n, mu, i):
                return RG.phi(T, n, mu, i) - delta_mu*n

            nmid = minmax_RG.maximize(newphi,T,nvapor,nliquid,mu,iterations)
            phi_mid = RG.phi(T, nmid, mu, iterations)
            print '      nmid,phi(nmid)',nmid,phi_mid

            nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,mu,iterations)
            nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,mu,iterations)
            print '      nvap,phi(nvap)',nvapor,phi_vapor
            print '      nliq,phi(nliq)',nliquid,phi_liquid
            print '\n'

        print '    left while loop'
        fout.write(str(T))
        fout.write('  ')
        fout.write(str(nvapor))
        fout.write('  ')
        fout.write(str(nmid))
        fout.write('  ')
        fout.write(str(nliquid))
        fout.write('  ')
        fout.write(str(phi_vapor))
        fout.write('  ')
        fout.write(str(phi_mid))
        fout.write('  ')
        fout.write(str(phi_liquid))
        fout.write('  ')
        fout.write(str(mu))
        fout.write('\n')
        sys.stdout.flush();
        print '   T, etaVap, etaLiq, mu',T,nvapor/(RG.sigma**3*np.pi/6),nliquid/(RG.sigma**3*np.pi/6),mu
        print '\n'
    fout.close()

if __name__ == '__main__':
    # Number of iterations desired
    # 0 iterations is TPT only, with no RG corrections
    num_iterations = 0

    # define the colors/symbols for plotting
    colors = np.array(['b-','g-','r-'])

    for i in range(num_iterations+1): # remember, range(x) only goes up to x-1
        print 'running npart_RG for iterations =',i

        # run npart
        npart(i)

        # Plot Coexistence
        # Read in data
        data = np.loadtxt('figs/coexistence_RG-i%d-out.dat'%i)

        T = data[:,0]
        etavapor = data[:,1]*np.pi*RG.sigma**3/6
        etaliquid = data[:,2]*np.pi*RG.sigma**3/6

    # Labels, etc.
    plt.xlabel(r'$\eta$')
    plt.ylabel('T')
    plt.legend(loc=0)
    plt.title('Liquid-Vapor Coexistence')

    # plt.savefig('figs/coexistence-RG-v2.pdf')
    plt.show()
