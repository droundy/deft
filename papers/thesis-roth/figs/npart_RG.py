from __future__ import division
import RG
import minmax_RG
import numpy as np
import time
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import forte

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

# input iteration depth as a parameter

def npart(iterations):
    # Initial conditions; dependent on you system

    T = 0.1
    nparticular = 0.08/(RG.sigma**3*np.pi/6) # using finite differences df_dn
    # nparticular = 0.155/(RG.sigma**3*np.pi/6) # using analytic df_nd

    # T = 0.68
    # nparticular = 0.0414/(RG.sigma**3*np.pi/6)

    N = 20 # For publication plots, make this bigger (40? 80? 100? You decide!)
    Tc = 1.33
    Tlow = T

    # Bounds of minimization.
    # Use range that works for many temperatures
    # Left (vapor)
    a_vap = 1e-10/(RG.sigma**3*np.pi/6)
    c_vap = nparticular
    # Right (liquid)
    a_liq = nparticular
    c_liq = 0.75/(RG.sigma**3*np.pi/6)
    ###############################

    # Open file for output
    fout = open('figs/npart_RG-i%d-out.dat'%iterations,'w')

    # label the columns of the output
    fout.write('#T     nvapor     nliquid       phi(nvap)        phi(nliq)         nparticular\n')#       phi_avg')
    sys.stdout.flush()

    # Do first temperature before the loop
    nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,nparticular,iterations)
    print 'initial nvap,phi_vap',nvapor,phi_vapor
    nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,nparticular,iterations)
    print 'initial nliq,phi_liq',nliquid,phi_liquid
    print 'initial npart',nparticular*(RG.sigma**3*np.pi/6)

    #while T < 1.4:
    for j in xrange(0,N+1):
        T = (Tc - Tlow)*(1 - ((N-j)/N)**4) + Tlow

        fout.flush()
        # Starting point for new nparticular is abscissa of max RG.phi with old nparticular
        nparticular = minmax_RG.maximize(RG.phi,T,nvapor, nliquid, nparticular,iterations)

        # I'm looking at the minima of RG.phi
        c_vap = nparticular
        a_liq = nparticular

        nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,c_vap,a_vap,nparticular,iterations)
        nliquid,phi_liquid = minmax_RG.minimize(RG.phi,T,a_liq,c_liq,nparticular,iterations)

        tol = 1e-5
        fnpart = RG.phi(T, nparticular, nparticular, iterations)

        # Compare the two minima in RG.phi
        while np.fabs(phi_vapor - phi_liquid)/np.fabs(fnpart) > tol:

            delta_mu = (phi_liquid - phi_vapor)/(nliquid - nvapor)

            def newphi(T, n, npart, i):
                return RG.phi(T, n, npart, i) - delta_mu*n

            nparticular = minmax_RG.maximize(newphi,T,nvapor,nliquid,nparticular,iterations)

            fnpart = RG.phi(T, nparticular, nparticular, iterations)

            nvapor,phi_vapor = minmax_RG.minimize(RG.phi,T,a_vap,c_vap,nparticular,iterations)
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
    fout.close()

if __name__ == '__main__':

    # I also want to know the time it takes
    timing = open('figs/RG-timing.dat','w')
    timing.write('# iterations time (s)\n')
    sys.stdout.flush()

    # Number of iterations desired
    num_iterations = 2

    # define the colors/symbols for plotting
    colors = np.array(['b-','g-','r-'])

    for i in range(num_iterations):
      print 'running npart_RG for iterations =',i
      timing.flush()
      # note current time
      time0 = time.clock()

      # run npart
      npart(i)

      # calculate elapsed time
      elapsed = time.clock() - time0

      # write to file
      timing.write(str(i)+' '+str(elapsed)+'\n')
      sys.stdout.flush()

      # Plot Coexistence
      # Read in data
      data = np.loadtxt('figs/npart_RG-i%d-out.dat'%i)

      T = data[:,0]
      etavapor = data[:,1]*np.pi*RG.sigma**3/6
      etaliquid = data[:,2]*np.pi*RG.sigma**3/6

      # Plot Forte's data
      plt.plot(etavapor, T, colors[i],label='RG '+r'$i=$ '+'%d'%i)
      plt.plot(etaliquid, T, colors[i])
      plt.plot(forte.eta,forte.T, label='Forte')

      # Labels, etc.
      plt.xlabel(r'$\eta$')
      plt.ylabel('T')
      plt.legend(loc=0)
      plt.title('Liquid-Vapor Coexistence')

      plt.savefig('figs/coexistance-RG.pdf')

    timing.close()
