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
# deft/papers/thesis-roth/figs
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

N = 100 # For publication plots, make this bigger (100? 500? 100000? You decide!)

# input iteration depth as a parameter
def npart(iterations):
    # Initial conditions; dependent on you system

    T = 0.5
    mu = -5.751 # I found this by plotting. It is good for T=0.5; i=0
    nmid = 0.0439651327587 # I found this from running the code with the above mu
    mu = -RG.df_dn(T, nmid, iterations)

    Tc = 1.3295 # Rough estimate based on previous plots
    Tlow = T

    # Bounds of minimization.
    # Use range that works for many temperatures
    # Outside bounds
    a_vap = 1e-10/(RG.sigma**3*np.pi/6)
    c_liq = 0.75/(RG.sigma**3*np.pi/6)
    # Use the central max as the inside bound
    nmid = 0.2/(RG.sigma**3*np.pi/6) # ad-hoc for now
    # nmid = minmax_RG.maximize(RG.phi,T,a_vap,c_liq,mu,iterations)

    ###############################

    # Open file for output
    fout = open('figs/coexistence_RG-i%d-out.dat'%iterations,'w')

    # label the columns of the output
    fout.write('#T     nvapor    nmid     nliquid       phi(nvap)       phi(nmid)         phi(nliq)         mu\n')
    sys.stdout.flush()

    # Do first temperature before the loop
    nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,nmid,mu,iterations)
    nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,nmid,c_liq,mu,iterations)
    print '  initial nvap,phi_vap',nvapor,phi_vapor
    print '  initial nmid,phi_mid',nmid,RG.phi(T,nmid,mu,iterations)
    print '  initial nliq,phi_liq',nliquid,phi_liquid
    print '  initial mu',mu

    #while T < 1.4:
    for j in xrange(0,N+1):
        # Temp points get closer as we near the critical point
        T = (Tc - Tlow)*(1 - ((N-j)/N)**4) + Tlow
        print '  T =',T

        fout.flush()

        # I'm looking at the minima of RG.phi
        # Use the local max between vapor and liquid to set the boundary for each minimization
        #nmid = minmax_RG.maximize(RG.phi,T,nvapor, nliquid, mu,iterations)
        mu = -RG.df_dn(T, nmid, iterations)
        #a_vap = max(nmid - 2.0*(nmid - nvapor), a_vap)
        #c_liq = nmid + 2.0*(nliquid - nmid)

        nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,nmid,mu,iterations)
        nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,nmid,c_liq,mu,iterations)
        phi_mid = RG.phi(T, nmid, mu, iterations)

        tol = 1e-9

        # Compare the two minima in RG.phi
        # print '    entering while loop'
        while np.fabs(phi_vapor - phi_liquid)/np.fabs(phi_mid) > tol: # np.fabs casts output as float
            # print '    whilecond =',np.fabs(phi_vapor - phi_liquid)/np.fabs(phi_mid)

            delta_mu = (phi_liquid - phi_vapor)/(nliquid - nvapor)
            # print '      delta_mu=',delta_mu

            # Change mu
            mu += delta_mu

            # find new values for nvap, nmid, nliq and phi_vap, phi_mid, phi_liq
            nmid, phi_mid = minmax_RG.maximize(RG.phi,T,nvapor,nliquid,mu,iterations), RG.phi(T, nmid, mu, iterations)
            nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,nmid,mu,iterations)
            nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,nmid,c_liq,mu,iterations)
            # print '      nvap,phi(nvap)',nvapor,phi_vapor
            # print '        etavap',nvapor*RG.sigma**3*np.pi/6
            # print '      nmid,phi(nmid)',nmid,phi_mid
            # print '        etamid',nmid*RG.sigma**3*np.pi/6
            # print '      nliq,phi(nliq)',nliquid,phi_liquid
            # print '        etaliq',nliquid*RG.sigma**3*np.pi/6
            # print '\n'

        if j > N - 20:
          plt.figure()
          delta_n = nliquid - nvapor
          nvals = np.arange(max(nvapor - delta_n, a_vap), min(nliquid + delta_n, c_liq), 0.01*delta_n)
          phivals = np.zeros_like(nvals)
          for i in range(len(nvals)):
            phivals[i] = RG.phi(T,nvals[i],mu,iterations)
          plt.plot(nvals*np.pi*RG.sigma**3/6, phivals, 'b-')
          plt.plot(nliquid*np.pi*RG.sigma**3/6, RG.phi(T,nliquid,mu,iterations), 'o')
          plt.plot(nvapor*np.pi*RG.sigma**3/6, RG.phi(T,nvapor,mu,iterations), 'o')
          plt.plot(nmid*np.pi*RG.sigma**3/6, RG.phi(T,nmid,mu,iterations), '+')
          plt.title('T = %.14g' % T)

        if nmid == nvapor or nmid == nliquid:
          print 'I have achieved silliness!'
          #break
        # print '    left while loop'
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
        # print '   T, etaVap, etaLiq, mu',T,nvapor/(RG.sigma**3*np.pi/6),nliquid/(RG.sigma**3*np.pi/6),mu
        # print '\n'
    fout.close()

if __name__ == '__main__':
    # Number of iterations desired
    # 0 iterations is TPT only, with no RG corrections
    num_iterations = 0

    # define the colors/symbols for plotting
    colors = np.array(['b-','g-','r-'])

    plt.figure(1)
    for i in range(num_iterations+1): # remember, range(x) only goes up to x-1
        print 'running npart_RG for iteration =',i

        # run npart
        npart(i)

        # Plot Coexistence
        # Read in data
        data = np.loadtxt('figs/coexistence_RG-i%d-out.dat'%i)

        T = data[:,0]
        etavapor = data[:,1]*np.pi*RG.sigma**3/6
        etaliquid = data[:,3]*np.pi*RG.sigma**3/6

        # Plot data
        plt.figure(1)
        plt.plot(etavapor, T, colors[i],label='RG '+r'$i=$ '+'%d'%i)
        plt.plot(etaliquid, T, colors[i])

    # Labels, etc.
    plt.xlabel(r'$\eta$')
    plt.ylabel('T')
    plt.legend(loc=0)
    plt.title('Liquid-Vapor Coexistence')

    # plt.savefig('figs/coexistence-RG-v2-80N.pdf')
    plt.savefig('figs/coexistence-RG-v2-%dN.pdf'%N)
    plt.show()
