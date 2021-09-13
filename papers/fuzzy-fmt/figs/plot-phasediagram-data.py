#!/usr/bin/python3

# This program produces temperature vs density, and pressure vs temperature
# phase diagrams as well as plots of P-vs-T at fixed n, P-vs_V at fixed T,
# P-vs-n at fixed T, and T-vs-n at fixed P from data produced by
# compute-phasediagram-data.py and stored in kT-%.3f.dat data files
# with the 3-column format: density, hfe, cfe.

# This program is used in conjunction compute-phasediagram-data.py which
# takes data from *best_tensor.dat data files generated by new-melting.cpp
# and creates the kT-%.3f.dat files used by this program.

# NOTE: Run this program from directory deft/papers/fuzzy-fmt
# with command: ./figs/plot-phasediagram-data.py  [show (to show plots)]

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
p_at_freezing = []
n_homogeneous_at_freezing = []
n_crystal_at_freezing = []
kT_homogeneous_at_freezing = []
kT_crystal_at_freezing = []
kT_in_plot = []

kT_data = []
density_data = []  # index corresponds to kT
pressure_data = []  # index corresponds to kT

our_kTs = (0.5, 1, 1.5, 2, 2.5, 3)

homogeneous_data = np.loadtxt('figs/homogeneous.dat')
homogeneous_temperature = homogeneous_data[0, 1:]
print(homogeneous_temperature)
homogeneous_density = homogeneous_data[1:, 0]
homogeneous_pressures = homogeneous_data[1:, 1:]

# for kT in (0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3):
# for kT in (0.6, 0.8, 1, 1.2, 1.4, 1.8, 2.2, 2.6, 3): #for paper
for kT in our_kTs:  # for paper
    # for kT in (0.6, 0.8, 1, 1.2, 1.4, 1.8, 2.2, 2.6, 3, 10, 14, 22, 30): #added higher temps
    n = []
    invn = []
    hfe = []
    cfe = []

    f = ('figs/crystal-data/kT-%.3f.dat' % kT)
    data = np.loadtxt(f)

    n = data[:, 0]  # density
    invn = 1/n
    hfe = data[:, 1]  # homogeneous free energy/atom
    cfe = data[:, 2]  # crystal free energy/atom

    hfe = np.array(hfe)
    cfe = np.array(cfe)
    invn = np.array(invn)
    n = np.array(n)

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
    #print('residuals', A[1])
    #print('coeff', coeff)
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
                # print(m2) #debug ASK!
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

    p_at_freezing.append(p_inter)
    n_homogeneous_at_freezing.append(1/invnh)
    n_crystal_at_freezing.append(1/invnc)

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

    kT_data.append(kT)  # holds all values of kT in a list
    density_data.append(actual_density)
    pressure_data.append(actual_pressure)

n_homogeneous_at_freezing = np.array(n_homogeneous_at_freezing)
n_crystal_at_freezing = np.array(n_crystal_at_freezing)
p_at_freezing = np.array(p_at_freezing)

plt.figure(figsize=myfigsize)
plt.fill_betweenx(kT_data, n_homogeneous_at_freezing,
                  n_crystal_at_freezing, color='#eeeeee')

# Plot T vs n  at constant P
for p in [2, 5, 10, 20]:
    n_mid_at_p_list = []
    kT_at_p_list = []
    for i in range(0, len(kT_data)):  # number of temperatures kT
        # number of elements of n at some kT
        for j in range(0, len(density_data[i])-1):
            if pressure_data[i][j] < p < pressure_data[i][j+1]:
                phi = pressure_data[i][j+1]
                plo = pressure_data[i][j]
                nhi = density_data[i][j+1]
                nlo = density_data[i][j]
                n_mid_at_p_list.append(
                    (nlo*(phi - p) + nhi*(p - plo))/(phi - plo))
                kT_at_p_list.append(kT_data[i])

    plt.plot(n_mid_at_p_list, kT_at_p_list, '.-', label='P=%g' % p)
plt.legend(loc='best')
plt.xlabel('$n^*$')
plt.ylabel('$T^*$')
plt.savefig('./figs/T-vs-n_at_fixed_P.pdf', transparent=True)

# - OR - uncomment the plot you want

# Plot n vs T  at constant P
#plt.plot(kT_at_p_list, n_mid_at_p_list, '.-', label= 'P=%g' % p)
#plt.title("Number Density vs Temperature at fixed Pressure")
# plt.legend(loc='best')
#plt.ylabel('Number Density')
# plt.xlabel('Temperature')

plt.figure(figsize=myfigsize)

plt.fill_betweenx(p_at_freezing, n_homogeneous_at_freezing,
                  n_crystal_at_freezing, color='#eeeeee')
for i in range(len(kT_data)):
    if kT_data[i] in [0.1, 0.2, 0.5, 1.0] or True:
        # Plot P vs n  at constant kT
        plt.plot(density_data[i], pressure_data[i],
                 label=r'$T^*=%g$' % kT_data[i])
plt.legend(loc='best')
#plt.ylim(0, 26)
#plt.ylim(0, 500)
plt.ylim(0, 45)
#plt.xlim(0, 1.1)
#plt.xlim(0, 1.8)
plt.xlim(0, 1.1)
plt.xlabel('$n^*$')
plt.ylabel('$p^*$')
plt.savefig('./figs/p-vs-n_at_fixed_T.pdf', transparent=True)

plt.figure('p vs volume', figsize=myfigsize)

temperatures_to_isotherm = [0.5, 1, 1.5, 2.0, 2.5, 3, 4, 5, 6, 8, 10, 15, 20]
temperatures_to_isotherm = [0.5, 1, 2.0, 5, 10, 15, 20, 30]

# Plot P vs 1/n (or V) at constant kT
plt.fill_betweenx(p_at_freezing, 1/n_homogeneous_at_freezing,
                  1/n_crystal_at_freezing, color='#eeeeee')
for i in range(len(kT_data)):
    if kT_data[i] in temperatures_to_isotherm:
        plt.plot(1/density_data[i], pressure_data[i],
                 color=styles.color_from_kT(kT_data[i]),
                 label='kT=%g' % kT_data[i])

for i in range(99, len(homogeneous_temperature)):
    if homogeneous_temperature[i] in [T for T in temperatures_to_isotherm if T not in kT_data]:
        print(f'temp[{i}] = {homogeneous_temperature[i]}')
        plt.plot(1/homogeneous_density, homogeneous_pressures[:, i], '--',
                 color=styles.color_from_kT(homogeneous_temperature[i]),
                 label='kT=%g' % homogeneous_temperature[i])


def R_BH(kT):
    eps = 1.0
    sigma = 1.0  # /* Let's define sigma == 1 for this one? */
    R = sigma*2.0**(1.0/6.0)
    N = 10000
    bh_diameter = 0
    dr = R/N
    beta = 1.0/kT
    for r_cur in np.arange(dr/2, R, dr):
        s6 = (sigma/r_cur)**6
        bh_diameter += (1 - np.exp(-beta*(4*eps*(s6*s6 - s6) + eps)))*dr
    return bh_diameter/2


def p_carnahan_starling(density, kT):
    eta = density*4*np.pi*R_BH(kT)**3/3
    return density*kT*(1 + eta + eta**2 - eta**3)/(1 - eta)**3


for kT in temperatures_to_isotherm:
    V_BH = 4*np.pi*R_BH(kT)**3/3
    V_freezing = V_BH/0.494  # This is the volume at which the fluid should freeze
    V_melting = V_BH/0.545  # This is the volume at which the solid should melt

    my_V = np.arange(V_melting, 3, 0.01)
    my_p = p_carnahan_starling(1/my_V, kT)
    my_p[my_V < V_freezing] = p_carnahan_starling(1/V_freezing, kT)

    plt.plot(my_V, my_p, ':',
             color=styles.color_from_kT(kT),
             # label='BH kT=%g' % kT
             )

# for temp in our_kTs:
#     volumes = []
#     pressures = []
#     for rd in np.arange(1.00, 0.25, -0.1):
#         fname = 'figs/mcfcc-%.4f-%.4f.dat.prs' % (rd, temp)
#         if os.path.exists(fname):
#             p = np.loadtxt(fname)
#             # we account below for the fact that the MC code does not spit
#             # out a reduced pressure, so we need to do a conversion.
#             p = p*2.0**(5/2.0)
#             pressures.append(p)
#             volumes.append(1/rd)
#     plt.plot(volumes, pressures, '.',
#              color=styles.color_from_kT(temp),
#              # label='MC %g' % temp
#              )

# for T in temperatures_to_isotherm:
#     V_zeno_T = []
#     p_zeno_T = []
#     for V in np.arange(0.9, 5, 0.01):
#         fname = 'data/mc/z-wca-108-%.2f-p-vs-T.dat' % V
#         if os.path.exists(fname):
#             data = np.loadtxt(fname)
#             p_zeno = data[:, 1]
#             T_zeno = data[:, 0]
#             for i in range(len(p_zeno)):
#                 if T_zeno[i] == T:
#                     p_zeno_T.append(p_zeno[i])
#                     V_zeno_T.append(V)
#     plt.plot(V_zeno_T, p_zeno_T, '.-.',
#                 color=styles.color_from_kT(T),
#                 linewidth=0.5,
#                 )

for T in temperatures_to_isotherm:
    V_zeno_T = []
    p_zeno_T = []
    for V in np.arange(0.9, 5, 0.01):
        fname = 'data/mc/z-wca-256-%.2f-p-vs-T.dat' % V
        if os.path.exists(fname):
            data = np.loadtxt(fname)
            p_zeno = data[:, 1]
            T_zeno = data[:, 0]
            for i in range(len(p_zeno)):
                if T_zeno[i] == T:
                    p_zeno_T.append(p_zeno[i])
                    V_zeno_T.append(V)
    plt.plot(V_zeno_T, p_zeno_T, '+-.',
                color=styles.color_from_kT(T),
                linewidth=0.5,
                )

# TODO Items week after May 19 2021
#
# - Change this code to make the colors the same at each temperature.  Show minimal legends.
#
# - Maybe try updating other visualizations to see if they are better
#
# - Try changing V range, p range, and T choices for clarity
#
# - Start looking again at the paper text.  Does anything need writing?
#
# - David: see about getting the new MC data into here?

# for temp in our_kTs:
#     volumes = []
#     pressures = []
#     for rd in np.arange(1.8, 0.25, -0.01):
#         fname = 'data/mc/ff-%.4g_temp-%.1f-press.dat' % (rd, temp)
#         if os.path.exists(fname):
#             p = np.loadtxt(fname)
#             # we account below for the fact that the MC code does not spit
#             # out a reduced pressure, so we need to do a conversion.
#             # ????
#             # p = p*2.0**(5/2.0)

#             # Is it p or p_excess?!
#             # p = p + rd*temp
#             pressures.append(p)
#             volumes.append(1/rd)
#     if len(volumes) > 0:
#         plt.plot(volumes, pressures, 'x',
#                  color=styles.color_from_kT(temp),
#                  # label='MC %g' % temp
#                  )


plt.plot([], [], 'k--', label='SFMT')
plt.plot([], [], 'k.-.', linewidth=0.5, label='Monte Carlo - Zeno')
# plt.plot([], [], 'kx', label='Monte Carlo - Chris')
plt.plot([], [], 'k:', label='Barker Henderson')
plt.legend(loc='upper right')
plt.ylim(0, 60)
plt.xlim(0.9, 2)
plt.xlabel('Volume per atom')
plt.ylabel('$p^*$')
plt.savefig('./figs/p-vs-V_at_fixed_T.pdf', transparent=True)

plt.figure(figsize=myfigsize)

densities_to_plot = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

# --------------NEW
# Plot P vs T  at constant n
# for n in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]:    #densities to show on the plot
# for n in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1]:  #densities to show on the plot
for n in densities_to_plot:  # densities to show on the plot in paper
    # for n in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:  #densities to show on the plot in thesis
    p_mid_at_n_list = []
    kT_at_n_list = []
    for i in range(0, len(kT_data)):  # number of temperatures kT
        # number of elements of P at some kT
        for j in range(0, len(pressure_data[i])-1):
            if density_data[i][j] < n < density_data[i][j+1]:
                phi = pressure_data[i][j+1]
                plo = pressure_data[i][j]
                nhi = density_data[i][j+1]
                nlo = density_data[i][j]
                p_mid_at_n_list.append(
                    (plo*(nhi - n) + phi*(n - nlo))/(nhi - nlo))
                kT_at_n_list.append(kT_data[i])

    plt.plot(kT_at_n_list, p_mid_at_n_list,
             styles.density_color(n)+'-', label='n=%g' % n)

for n in densities_to_plot:
    fname = 'data/mc/z-wca-108-n%.2f-p-vs-T.dat' % n
    if os.path.exists(fname):
        data = np.loadtxt(fname)
        p_zeno = data[:, 1]
        T_zeno = data[:, 0]
        plt.plot(T_zeno, p_zeno, '-.',
                    color=styles.density_color(n),
                    linewidth=0.5)
for n in densities_to_plot:
    fname = 'data/mc/z-wca-256-n%.2f-p-vs-T.dat' % n
    if os.path.exists(fname):
        data = np.loadtxt(fname)
        p_zeno = data[:, 1]
        T_zeno = data[:, 0]
        plt.plot(T_zeno, p_zeno, '.--',
                    color=styles.density_color(n),
                    linewidth=0.5)

for rd in np.arange(1.00, 0.25, -0.1):
    temps = []
    pressures = []
    # print 'could use:', glob.glob('figs/mcfcc-%.4f-*.dat.prs' % (rd))
    for temp in np.arange(0.05, 3.05, 0.05):
        fname = 'figs/mcfcc-%.4f-%.4f.dat.prs' % (rd, temp)
        if os.path.exists(fname):
            temps.append(temp)
            p = np.loadtxt(fname)
            # we account below for the fact that the MC code does not spit
            # out a reduced pressure, so we need to do a conversion.
            pressures.append(p*2.0**(5/2.0))
    if len(temps) > 0:
        plt.plot(temps, pressures, styles.density_color(rd)+'--')
        # if len(np.argwhere(np.abs(Nlabel - rd) < 0.001)) == 1:
        #   print('We have MC for rd', rd)
        #   pressures = np.array(pressures)
        #   temps = np.array(temps)
        # i = np.argwhere(np.abs(Nlabel - rd) < 0.001)[0][0]
        # j = np.argwhere(temps > Tlabel[i])[0][0]
        # MCPlabel[i] = pressures[j]

bhdata = np.loadtxt('figs/homogeneous-bh.dat')
bhtemp = bhdata[0, 1:]
bhnred = bhdata[1:, 0]
bhpressures = bhdata[1:, 1:]
bhnans = bhpressures != bhpressures
bhpressures[bhnans] = 1e30
bhtemperatures, bhn_reduced = np.meshgrid(bhtemp, bhnred)

# Plot Barker-Henderson results
for i in range(9, len(bhn_reduced[:, 0]), 10)[:10]:
    plt.plot(bhtemperatures[i, :], bhpressures[i, :],
             styles.density_color(bhn_reduced[i, 0])+':')

plt.legend(loc='best')
plt.ylabel('$p^*$')
plt.xlabel('$T^*$')
plt.xlim(0, 3)  # for paper
# plt.xlim(0.5,3) #for thesis
plt.ylim(0, 35)
plt.tight_layout()
plt.savefig('./figs/p-vs-T_at_fixed_n.pdf', transparent=True)

# - OR - uncomment the plot you want
# Plot T vs P  at constant n
#plt.plot(kT_at_n_list, p_mid_at_n_list, '.-', label= 'n=%g' % n)
#plt.title("Temperature vs Pressure at fixed n")
# plt.legend(loc='best')
# plt.xlabel('Pressure')
# plt.ylabel('Temperature')
# --------------end NEW

plt.figure(figsize=myfigsize)

# Temperature vs Density Phase Diagram
plt.plot(n_homogeneous_at_freezing, kT_data, label='fluid', color='red')
plt.plot(n_crystal_at_freezing, kT_data, label='solid', color='blue')
plt.fill_betweenx(kT_data, .2, n_homogeneous_at_freezing, color='red')
plt.fill_betweenx(kT_data, n_homogeneous_at_freezing,
                  n_crystal_at_freezing, color='gray')
#plt.fill_betweenx(kT_data, n_crystal_at_freezing, 1.6, color='blue')
plt.fill_betweenx(kT_data, n_crystal_at_freezing, 1.8, color='blue')
# plt.legend(loc='best')
plt.xlabel('$n^*$')
plt.ylabel('$T^*$')

##plt.plot([0.88, 0.90, 0.91, 0.92, 1.04, 1.12],[0.7, 0.8, 0.9, 1.0, 2.0, 3.0], label='chris_l', color='green')
##plt.plot([0.96, 0.98, 0.99, 1.00, 1.11, 1.19],[0.7, 0.8, 0.9, 1.0, 2.0, 3.0], label='chris_s', color='green')
#plt.plot([0.88, 0.90, 0.91, 0.92, 1.04, 1.12, 1.24, 1.44],[0.7, 0.8, 0.9, 1.0, 2.0, 3,5,10], label='chris_l', color='green')
#plt.plot([0.96, 0.98, 0.99, 1.00, 1.11, 1.19, 1.31, 1.51],[0.7, 0.8, 0.9, 1.0, 2.0, 3, 5, 10], label='chris_s', color='green')
plt.plot([0.88, 0.90, 0.91, 0.92, 1.04, 1.12], [0.7, 0.8, 0.9,
         1.0, 2.0, 3], '-', label='$MC_f$', color='white')
plt.plot([0.96, 0.98, 0.99, 1.00, 1.11, 1.19], [0.7, 0.8, 0.9,
         1.0, 2.0, 3], '-', label='$MC_s$', color='white')
plt.legend()
plt.xlim(0.2, 1.8)
plt.ylim(0.5, 3)

plt.tight_layout()
plt.savefig('./figs/Phase_Diagram_of_T_vs_n.pdf', transparent=True)

plt.figure(figsize=myfigsize)

# Pressure vs Temperature Phase Diagram
plt.fill_between(kT_data, 0, p_at_freezing, label='fluid', color='red')
# plt.fill_between(kT_data, p_at_freezing, 50, color='blue')    #FIX - change 30
#plt.fill_between(kT_data, p_at_freezing, 1500, color='blue')
plt.fill_between(kT_data, p_at_freezing, 80, label='solid', color='blue')
plt.plot(kT_data, p_at_freezing, color='black')
#plt.ylim(0, 40)
# plt.xlim(kT_data.min(), kT_data.max())     #FIX!
plt.xlabel('$T^*$')
plt.ylabel('$p^*$')
#plt.plot([0.7, 0.8,0.9,1.0,2.0,3.0], [6.24, 7.62, 8.78, 9.99, 25.5,43.8], label='chris_l', color='green')
##plt.plot([0.7, 0.8,0.9,1.0,2.0, 3, 5, 10], [6.24, 7.62, 8.78, 9.99, 25.5,43.8, 85.6, 210], label='chris_l', color='green')
plt.plot([0.7, 0.8, 0.9, 1.0, 2.0, 3], [6.24, 7.62, 8.78,
         9.99, 25.5, 43.8], label='MC', color='white')
plt.legend()
plt.ylim(0, 80)
plt.xlim(0.5, 3)

plt.tight_layout()
plt.savefig('./figs/Phase_Diagram_of_P_vs_T.pdf', transparent=True)

if 'show' in sys.argv:
    plt.show()
