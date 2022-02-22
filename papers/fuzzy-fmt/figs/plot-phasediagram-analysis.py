#!/usr/bin/python3

# This program produces plots of  P-vs-n at fixed T and P-vs_V at fixed T
# from data produced by compute-phasediagram-data.py 
# and stored in kT-%.3f-seed%g.dat data files
# with the 5-column format: density, hfe, cfe, hpressure-nm, cpressure-nm

# NOTE: Run this program from directory deft/papers/fuzzy-fmt
# with command: ./figs/plot-phasediagram-analysis.py  [show (to show plots)]

from __future__ import print_function, division

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import sys
import styles

from cycler import cycler

matplotlib.rc('text', usetex=True)
figscale = 1.5
myfigsize = (4*figscale, 3*figscale)

# pressure at freezing (at intersection of homogeneous and crystal plots)
p_at_freezing_seed1 = []
p_at_freezing_seed2 = []
p_at_freezing_seed3 = []
p_at_freezing_seed4 = []
p_at_freezing_seed5 = []
p_at_freezing_seed6 = []   #not really seed=6, rather, seed=1 at higher accuracy!
n_homogeneous_at_freezing_seed1 = []
n_homogeneous_at_freezing_seed2 = []
n_homogeneous_at_freezing_seed3 = []
n_homogeneous_at_freezing_seed4 = []
n_homogeneous_at_freezing_seed5 = []
n_homogeneous_at_freezing_seed6 = []  #not really seed=6, rather, seed=1 at higher accuracy!
n_crystal_at_freezing_seed1 = []
n_crystal_at_freezing_seed2 = []
n_crystal_at_freezing_seed3 = []
n_crystal_at_freezing_seed4 = []
n_crystal_at_freezing_seed5 = []
n_crystal_at_freezing_seed6 = [] #not really seed=6, rather, seed=1 at higher accuracy!

kT_data_seed1 = []
kT_data_seed2 = []
kT_data_seed3 = []
kT_data_seed4 = []
kT_data_seed5 = []
kT_data_seed6 = []  #not really seed=6, rather, seed=1 at higher accuracy!
density_data_seed1 = []  # index corresponds to kT
density_data_seed2 = []  # index corresponds to kT
density_data_seed3 = []  # index corresponds to kT
density_data_seed4 = []  # index corresponds to kT
density_data_seed5 = []  # index corresponds to kT
density_data_seed6 = []  # index corresponds to kT  #not really seed=6, rather, seed=1 at higher accuracy!
pressure_data_seed1 = []  # index corresponds to kT
pressure_data_seed2 = []  # index corresponds to kT
pressure_data_seed3 = []  # index corresponds to kT
pressure_data_seed4 = []  # index corresponds to kT
pressure_data_seed5 = []  # index corresponds to kT
pressure_data_seed6 = []  # index corresponds to kT  #not really seed=6, rather, seed=1 at higher accuracy!
density_data_nm_seed1 = []  # index corresponds to kT
density_data_nm_seed2 = []  # index corresponds to kT
density_data_nm_seed3 = []  # index corresponds to kT
density_data_nm_seed4 = []  # index corresponds to kT
density_data_nm_seed5 = []  # index corresponds to kT
density_data_nm_seed6 = []  # index corresponds to kT  #not really seed=6, rather, seed=1 at higher accuracy!
pressure_data_nm_seed1 = []  # index corresponds to kT
pressure_data_nm_seed2 = []  # index corresponds to kT
pressure_data_nm_seed3 = []  # index corresponds to kT
pressure_data_nm_seed4 = []  # index corresponds to kT
pressure_data_nm_seed5 = []  # index corresponds to kT
pressure_data_nm_seed6 = []  # index corresponds to kT  #not really seed=6, rather, seed=1 at higher accuracy!

our_kTs = (0.5,3) #debug
seeds = (1, 2, 3, 4, 5, 6)

homogeneous_data = np.loadtxt('figs/homogeneous.dat')
homogeneous_temperature = homogeneous_data[0, 1:]
print(homogeneous_temperature)
homogeneous_density = homogeneous_data[1:, 0]
homogeneous_pressures = homogeneous_data[1:, 1:]

for seed in seeds:
  for kT in our_kTs:
    n = []
    invn = []
    hfe = []
    cfe = []
    hpressure_nm = []
    cpressure_nm = []
    actual_pressure_nm = []
    actual_density_nm = []

    f = ('figs/crystal-data/kT-%.3f-seed%g.dat' % (kT, seed))    #debug
    data = np.loadtxt(f)

    n = data[:, 0]  # density
    invn = 1/n
    hfe = data[:, 1]  # homogeneous free energy/atom
    cfe = data[:, 2]  # crystal free energy/atom
    hpressure_nm = data[:,3]
    cpressure_nm = data[:,4]

    hfe = np.array(hfe)
    cfe = np.array(cfe)
    invn = np.array(invn)
    n = np.array(n)
    
    for i in range(len(n)):
       if cfe[i] < hfe[i]:
           pressure_nm=cpressure_nm[i]
       if hfe[i] < cfe[i]:
           pressure_nm=hpressure_nm[i] 
       if hfe[i] == 0:     #debug only - delete when working
           pressure_nm=0  #debug only - delete when working            
       actual_pressure_nm.append(pressure_nm)
       actual_density_nm.append(n[i])  
                
    actual_pressure_nm = np.array(actual_pressure_nm) 
    actual_density_nm = np.array(actual_density_nm)               

    functions = np.vstack((np.ones_like(invn),
                           invn**-1,
                           invn**-2,
                           invn**-3,
                           invn**-4,
                           invn**-5,
                           invn**-6)).T
    pressure_functions = np.vstack((np.zeros_like(invn),
                                    invn**-2,
                                    2*invn**-3,
                                    3*invn**-4,
                                    4*invn**-5,
                                    5*invn**-6,
                                    6*invn**-7)).T
    A = np.linalg.lstsq(functions, cfe)
    coeff = A[0]
    fit_cfe = np.dot(functions, coeff)

    # pressure_skip is the number of points to skip over when
    # taking finite difference derivatives.
    pressure_skip = 1
    dhfe = hfe[pressure_skip:] - hfe[:-pressure_skip]
    dcfe = cfe[pressure_skip:] - cfe[:-pressure_skip]
    dinvn = invn[pressure_skip:] - invn[:-pressure_skip]
    mid_invn = invn[:-pressure_skip] + dinvn/2
    mid_hfe = hfe[:-pressure_skip] + dhfe/2
    mid_cfe = cfe[:-pressure_skip] + dcfe/2

    hpressure = -(dhfe/dinvn)  # for fixed N and Te
    cpressure = -(dcfe/dinvn)  # for fixed N and Te

    fit_p = np.dot(pressure_functions, coeff)

    mid_h_gibbs = mid_hfe + mid_invn*hpressure
    mid_c_gibbs = mid_cfe + mid_invn*cpressure
    fit_c_gibbs = fit_cfe + invn*fit_p

    # Find pressure at point of intersection
    def find_first_intersection(p1, g1, p2, g2):
        for i in range(1, len(g1)-1):
            m1 = (g1[i+1]-g1[i])/(p1[i+1]-p1[i])
            for j in range(1, len(g2)-1):
                m2 = (g2[j+1]-g2[j])/(p2[j+1]-p2[j])
                if m1 != m2:
                    P_inter = (g2[j] - m2*p2[j] - g1[i] + m1*p1[i])/(m1-m2)
                    if p1[i] < P_inter < p1[i+1] and p2[j] < P_inter < p2[j+1]:
                        g_inter = m1*P_inter+g1[i]-m1*p1[i]
                        if g1[i] < g_inter < g1[i+1] and g2[j] < g_inter < g2[j+1]:
                            return P_inter, g_inter

    p_inter, g_inter = find_first_intersection(
        hpressure, mid_h_gibbs, cpressure, mid_c_gibbs)
    pf_inter, gf_inter = find_first_intersection(
        hpressure, mid_h_gibbs, fit_p, fit_c_gibbs)

    # Find homogeneous and crystal densities at p_inter

    def find_densities(p_inter, pressure, invn):
        for i in range(1, len(pressure)-1):
            if pressure[i] > p_inter:
                pressureabove = pressure[i]
                invnabove = invn[i]
                pressurebelow = pressure[i-1]
                invnbelow = invn[i-1]
                m = (pressureabove-pressurebelow)/(invnabove-invnbelow)
                invn_inter = invnabove-((pressureabove-p_inter)/m)
                return invn_inter
    invnh = find_densities(p_inter, hpressure, mid_invn)
    invnc = find_densities(p_inter, cpressure, mid_invn)
    
    if seed ==1:
      p_at_freezing_seed1.append(p_inter)
      n_homogeneous_at_freezing_seed1.append(1/invnh)
      n_crystal_at_freezing_seed1.append(1/invnc)
    if seed ==2:
      p_at_freezing_seed2.append(p_inter)
      n_homogeneous_at_freezing_seed2.append(1/invnh)
      n_crystal_at_freezing_seed2.append(1/invnc)
    if seed ==3:
      p_at_freezing_seed3.append(p_inter)
      n_homogeneous_at_freezing_seed3.append(1/invnh)
      n_crystal_at_freezing_seed3.append(1/invnc)
    if seed ==4:
      p_at_freezing_seed4.append(p_inter)
      n_homogeneous_at_freezing_seed4.append(1/invnh)
      n_crystal_at_freezing_seed4.append(1/invnc)
    if seed ==5:
      p_at_freezing_seed5.append(p_inter)
      n_homogeneous_at_freezing_seed5.append(1/invnh)
      n_crystal_at_freezing_seed5.append(1/invnc)
    if seed ==6:
      p_at_freezing_seed6.append(p_inter)
      n_homogeneous_at_freezing_seed6.append(1/invnh)
      n_crystal_at_freezing_seed6.append(1/invnc)

    # compute the actual physical pressure as a function of density, and skip over coexistence
    actual_pressure = []
    actual_density = []
    for i in range(len(mid_invn)):
        if hpressure[i] >= p_inter:
            break  # if the pressure is too high, then we should just stop, since we have left the fluid
        actual_pressure.append(hpressure[i])
        actual_density.append(1/mid_invn[i])
    actual_pressure.append(p_inter)
    actual_density.append(1/invnh)
    actual_pressure.append(p_inter)
    actual_density.append(1/invnc)
    for i in range(len(mid_invn)):
        if cpressure[i] < 0 and mid_invn[i] <= invnc:
            # when the pressure is negative, we know we are in the crazy part where our dft fails.
            break
        if cpressure[i] > p_inter:
            actual_pressure.append(cpressure[i])
            actual_density.append(1/mid_invn[i])        
    actual_pressure = np.array(actual_pressure)
    actual_density = np.array(actual_density)
    
    # print (kT, p_inter, 1/invnh, 1/invnc)   #Use >> phase_diagram_data.dat (or phase_diagram_data-tensor.dat) to store data for reference

    if seed ==1:
       kT_data_seed1.append(kT)  # holds all values of kT in a list
       density_data_seed1.append(actual_density)
       pressure_data_seed1.append(actual_pressure)
       density_data_nm_seed1.append(actual_density_nm)
       pressure_data_nm_seed1.append(actual_pressure_nm)
    if seed ==2:
       kT_data_seed2.append(kT)  # holds all values of kT in a list
       density_data_seed2.append(actual_density)
       pressure_data_seed2.append(actual_pressure)
       density_data_nm_seed2.append(actual_density_nm)
       pressure_data_nm_seed2.append(actual_pressure_nm)
    if seed ==3:
       kT_data_seed3.append(kT)  # holds all values of kT in a list
       density_data_seed3.append(actual_density)
       pressure_data_seed3.append(actual_pressure)
       density_data_nm_seed3.append(actual_density_nm)
       pressure_data_nm_seed3.append(actual_pressure_nm) 
    if seed ==4:
       kT_data_seed4.append(kT)  # holds all values of kT in a list
       density_data_seed4.append(actual_density)
       pressure_data_seed4.append(actual_pressure)
       density_data_nm_seed4.append(actual_density_nm)
       pressure_data_nm_seed4.append(actual_pressure_nm)       
    if seed ==5:
       kT_data_seed5.append(kT)  # holds all values of kT in a list
       density_data_seed5.append(actual_density)
       pressure_data_seed5.append(actual_pressure)
       density_data_nm_seed5.append(actual_density_nm)
       pressure_data_nm_seed5.append(actual_pressure_nm) 
    if seed ==6:
       kT_data_seed6.append(kT)  # holds all values of kT in a list
       density_data_seed6.append(actual_density)
       pressure_data_seed6.append(actual_pressure)
       density_data_nm_seed6.append(actual_density_nm)
       pressure_data_nm_seed6.append(actual_pressure_nm) 
             
  if seed ==1:      
   n_homogeneous_at_freezing_seed1 = np.array(n_homogeneous_at_freezing_seed1)
   n_crystal_at_freezing_seed1 = np.array(n_crystal_at_freezing_seed1)
   p_at_freezing_seed1 = np.array(p_at_freezing_seed1)
  if seed ==2:      
   n_homogeneous_at_freezing_seed2 = np.array(n_homogeneous_at_freezing_seed2)
   n_crystal_at_freezing_seed2 = np.array(n_crystal_at_freezing_seed2)
   p_at_freezing_seed2 = np.array(p_at_freezing_seed2)
  if seed ==3:      
   n_homogeneous_at_freezing_seed3 = np.array(n_homogeneous_at_freezing_seed3)
   n_crystal_at_freezing_seed3 = np.array(n_crystal_at_freezing_seed3)
   p_at_freezing_seed3 = np.array(p_at_freezing_seed3)
  if seed ==4:      
   n_homogeneous_at_freezing_seed4 = np.array(n_homogeneous_at_freezing_seed4)
   n_crystal_at_freezing_seed4 = np.array(n_crystal_at_freezing_seed4)
   p_at_freezing_seed4 = np.array(p_at_freezing_seed4)
  if seed ==5:      
   n_homogeneous_at_freezing_seed5 = np.array(n_homogeneous_at_freezing_seed5)
   n_crystal_at_freezing_seed5 = np.array(n_crystal_at_freezing_seed5)
   p_at_freezing_seed5 = np.array(p_at_freezing_seed5)
  # if seed ==6:
   # n_homogeneous_at_freezing_seed6 = np.array(n_homogeneous_at_freezing_seed6)
   # n_crystal_at_freezing_seed6 = np.array(n_crystal_at_freezing_seed6)
   # p_at_freezing_seed6 = np.array(p_at_freezing_seed6)

# Plot P vs n  at constant kT
plt.figure('p vs n', figsize=myfigsize)
plt.fill_betweenx(p_at_freezing_seed1, n_homogeneous_at_freezing_seed1,
                  n_crystal_at_freezing_seed1, color='#eeeeee')
i=0  #kT=0.5
for seed in seeds:
 #print pressure data for kT=0.5 from differential for all seeds:
 if seed ==1:
            plt.plot(density_data_seed1[i], pressure_data_seed1[i], color='red',
                 label=r'$T^*=%g$~~seed%g' % (kT_data_seed1[i], seed))
 if seed ==2:
            plt.plot(density_data_seed2[i], pressure_data_seed2[i], color='green',
                 label=r'$T^*=%g$~~seed%g' % (kT_data_seed2[i], seed))
 if seed ==3:
            plt.plot(density_data_seed3[i], pressure_data_seed3[i], color='blue',
                 label=r'$T^*=%g$~~seed%g' % (kT_data_seed3[i], seed))
 if seed ==4:
            plt.plot(density_data_seed4[i], pressure_data_seed4[i], color='orange',
                 label=r'$T^*=%g$~~seed%g' % (kT_data_seed4[i], seed))
 if seed ==5:
            plt.plot(density_data_seed5[i], pressure_data_seed5[i], color='brown',
                 label=r'$T^*=%g$~~seed%g' % (kT_data_seed5[i], seed))
 # if seed ==6:  #not really seed=6, rather, seed=1 at higher accuracy!
            # plt.plot(density_data_seed6[i], pressure_data_seed6[i], color='black',
                 # label=r'$T^*=%g$~~seed1 at mcerr4' % (kT_data_seed6[i]))
                 
 # #print pressure data for kT=0.5 from new-melting for all seeds:
 # if seed ==1:
            # plt.plot(density_data_nm_seed1[i], pressure_data_nm_seed1[i], '--', color='red',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed1[i], seed))
 # if seed ==2:
            # plt.plot(density_data_nm_seed2[i], pressure_data_nm_seed2[i], '--', color='green',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed2[i], seed))
 # if seed ==3:
            # plt.plot(density_data_nm_seed3[i], pressure_data_nm_seed3[i], '--', color='blue',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed3[i], seed))
 # if seed ==4:
            # plt.plot(density_data_nm_seed4[i], pressure_data_nm_seed4[i], '--', color='orange',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed4[i], seed))
 # if seed ==5:
            # plt.plot(density_data_nm_seed5[i], pressure_data_nm_seed5[i], '--', color='brown',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed5[i], seed))
 # # if seed ==6:   #not really seed=6, rather, seed=1 at higher accuracy!
            # # plt.plot(density_data_nm_seed6[i], pressure_data_nm_seed6[i], '--', color='black',
                 # # label=r'$T^*=%gnm$~~seed1 at mcerr4' % (kT_data_seed6[i]))

# i=1  #kT=3
# for seed in seeds:
 # #print pressure data for kT=3 from differential for all seeds:
 # if seed ==1:
            # plt.plot(density_data_seed1[i], pressure_data_seed1[i], color='purple',
                 # label=r'$T^*=%g$~~seed%g' % (kT_data_seed1[i], seed))
 # if seed ==2:
            # plt.plot(density_data_seed2[i], pressure_data_seed2[i], color='xkcd:forest green',
                 # label=r'$T^*=%g$~~seed%g' % (kT_data_seed2[i], seed))
 # if seed ==3:
            # plt.plot(density_data_seed3[i], pressure_data_seed3[i], color='xkcd:dark blue',
                 # label=r'$T^*=%g$~~seed%g' % (kT_data_seed3[i], seed))
 # if seed ==4:
            # plt.plot(density_data_seed4[i], pressure_data_seed4[i], color='xkcd:hot pink',
                 # label=r'$T^*=%g$~~seed%g' % (kT_data_seed4[i], seed))
 # if seed ==5:
            # plt.plot(density_data_seed5[i], pressure_data_seed5[i], color='xkcd:magenta',
                 # label=r'$T^*=%g$~~seed%g' % (kT_data_seed5[i], seed))
 # # if seed ==6:   #not really seed=6, rather, seed=1 at higher accuracy!
            # # plt.plot(density_data_seed6[i], pressure_data_seed6[i], color='black',
                 # # label=r'$T^*=%g$~~seed1 at mcerr4' % (kT_data_seed6[i]))
                 

 # #print pressure data for kT=3 from new-melting for all seeds:
 # if seed ==1:
            # plt.plot(density_data_nm_seed1[i], pressure_data_nm_seed1[i], '--', color='purple',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed1[i], seed))
 # if seed ==2:
            # plt.plot(density_data_nm_seed2[i], pressure_data_nm_seed2[i], '--', color='xkcd:forest green',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed2[i], seed))
 # if seed ==3:
            # plt.plot(density_data_nm_seed3[i], pressure_data_nm_seed3[i], '--', color='xkcd:dark blue',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed3[i], seed))
 # if seed ==4:
            # plt.plot(density_data_nm_seed4[i], pressure_data_nm_seed4[i], '--', color='xkcd:hot pink',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed4[i], seed))
 # if seed ==5:
            # plt.plot(density_data_nm_seed5[i], pressure_data_nm_seed5[i], '--', color='xkcd:magenta',
                 # label=r'$T^*=%gnm$~~seed%g' % (kT_data_seed5[i], seed))
 # # if seed ==6:   #not really seed=6, rather, seed=1 at higher accuracy!
            # # plt.plot(density_data_nm_seed6[i], pressure_data_nm_seed6[i], '--', color='black',
                 # # label=r'$T^*=%gnm$~~seed1 at mcerr4' % (kT_data_seed6[i]))
    
plt.legend(loc='best')
plt.ylim(0, 45)
plt.xlim(0, 1.1)
plt.xlabel('$n^*$')
plt.ylabel('$p^*$')
plt.savefig('./figs/p-vs-n_at_fixed_T.pdf', transparent=True)



# # Plot P vs 1/n (or V) at constant kT
# plt.figure('p vs volume', figsize=myfigsize)

# temperatures_to_isotherm = [0.5, 3]

# plt.fill_betweenx(p_at_freezing_seed1, 1/n_homogeneous_at_freezing_seed1, 
                  # 1/n_crystal_at_freezing_seed1, color='#eeeeee')
# for i in range(len(kT_data_seed1)):
    # if kT_data_seed1[i] in temperatures_to_isotherm:
        # plt.plot(1/density_data_seed1[i], pressure_data_seed1[i],
                 # color=styles.color_from_kT(kT_data_seed1[i]),
                 # label='kT=%g~~seed%g' % (kT_data_seed1[i], seed))

# for i in range(99, len(homogeneous_temperature)):
    # if homogeneous_temperature[i] in [T for T in temperatures_to_isotherm if T not in kT_data_seed1]:
        # print(f'temp[{i}] = {homogeneous_temperature[i]}')
        
        # plt.plot(1/homogeneous_density, homogeneous_pressures[:, i], '-',
                 # color=styles.color_from_kT(homogeneous_temperature[i]),
                 # label='kT=%g' % homogeneous_temperature[i])

# def R_BH(kT):
    # eps = 1.0
    # sigma = 1.0  # /* Let's define sigma == 1 for this one? */
    # R = sigma*2.0**(1.0/6.0)
    # N = 10000
    # bh_diameter = 0
    # dr = R/N
    # beta = 1.0/kT
    # for r_cur in np.arange(dr/2, R, dr):
        # s6 = (sigma/r_cur)**6
        # bh_diameter += (1 - np.exp(-beta*(4*eps*(s6*s6 - s6) + eps)))*dr
    # return bh_diameter/2


# def p_carnahan_starling(density, kT):
    # eta = density*4*np.pi*R_BH(kT)**3/3
    # return density*kT*(1 + eta + eta**2 - eta**3)/(1 - eta)**3


# for kT in temperatures_to_isotherm:
    # V_BH = 4*np.pi*R_BH(kT)**3/3
    # V_freezing = V_BH/0.494  # This is the volume at which the fluid should freeze
    # V_melting = V_BH/0.545  # This is the volume at which the solid should melt

    # my_V = np.arange(V_melting, 3, 0.01)
    # my_p = p_carnahan_starling(1/my_V, kT)
    # my_p[my_V < V_freezing] = p_carnahan_starling(1/V_freezing, kT)

    # plt.plot(my_V, my_p, ':',
             # color=styles.color_from_kT(kT),
             # # label='BH kT=%g' % kT
             # )

# for T in temperatures_to_isotherm:
    # V_zeno_T = []
    # p_zeno_T = []
    # for V in np.arange(0.9, 5, 0.01):
        # fname = 'data/mc/z-wca-256-%.2f-p-vs-T.dat' % V
        # if os.path.exists(fname):
            # data = np.loadtxt(fname)
            # p_zeno = data[:, 1]
            # T_zeno = data[:, 0]
            # for i in range(len(p_zeno)):
                # if T_zeno[i] == T:
                    # p_zeno_T.append(p_zeno[i])
                    # V_zeno_T.append(V)
    # plt.plot(V_zeno_T, p_zeno_T, '.-.',
                # color=styles.color_from_kT(T),
                # linewidth=0.5,
                # )
                
# plt.plot([], [], 'k-', label='SFMT')
# plt.plot([], [], 'k.-.', linewidth=0.5, label='Monte Carlo - Zeno')
# # plt.plot([], [], 'kx', label='Monte Carlo - Chris')
# plt.plot([], [], 'k:', label='Barker Henderson')
# plt.legend(loc='upper right')
# plt.ylim(0, 60)
# plt.xlim(0.9, 2)
# plt.xlabel('Volume per atom')
# plt.ylabel('$p^*$')
# plt.savefig('./figs/p-vs-V_at_fixed_T.pdf', transparent=True)

if 'show' in sys.argv:
    plt.show()
