#!/usr/bin/python3

#This program creates plots of: Free Energy/atom vs 1/Reduced Density
#                               Gibbs Free Energy/atom vs Pressure
#                               Pressure vs 1/Reduced Density
#for a given temperature from data in kT*best.dat (or kT*best_tensor.dat) files
#which are generated as output data files by figs/new-melting.cpp
#This program is the same plot-pressure.py prgram modified for plots labeled nicely for thesis 

#NOTE: Run this plot script from directory deft/papers/thesis-kirstie 
#with comand figs/plot-pressure-thesis.py --kT [temp] [OPTIONAL: --tensor  --showplots --saveplots]

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
#import argparse

#parser = argparse.ArgumentParser(description='Creates plots of: FE/atom vs 1/n, GibbsFE/atom vs P, and P vs 1/n')

#parser.add_argument('--kT', metavar='temperature', type=float,
#                    help='reduced temperature - REQUIRED')

#parser.add_argument('directory', metavar='directory', type=str,
#                    help='directory with data to plot') 

#parser.add_argument('--tensor', action='store_true',
#                    help='use data with tensor weight')

#parser.add_argument('--showplots', action='store_true',
#                    help='show plots')
                    
#parser.add_argument('--saveplots', action='store_true',
#                    help='will generate 3 plot .pdf files')

#args=parser.parse_args()

#kT=args.kT
kT=10

n = []
invn = []
hfe = []
cfe = []

#files = sorted(list(glob.glob('../fuzzy-fmt/newdata_tensor/phase-diagram4/kT%.3f_n*_best_tensor.dat' % kT)))  #works with fac, but data files must be in git
files = sorted(list(glob.glob('figs/plot-pressure-data/kT%.3f_n*_best_tensor.dat' % kT))) #works with fac, but data files must be in git

for f in files:
    data = np.loadtxt(f)
    n.append(data[1])     #density
    invn.append(1/data[1])
    hfe.append(data[4])   #homogeneous free energy/atom
    cfe.append(data[5])   #crystal free energy/atom
hfe = np.array(hfe)
cfe = np.array(cfe)
invn = np.array(invn)
n = np.array(n)

print(n)

# Generate data for Figure 1
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


dhfe=np.diff(hfe)  #Caution: depends on order of data files!
dcfe=np.diff(cfe)  #Caution: depends on order of data files!
dinvn=np.diff(invn)  #Caution: depends on order of data files!
mid_invn=invn[0:len(invn)-1]+dinvn/2
hpressure = -(dhfe/dinvn) #for fixed N and Te   
cpressure = -(dcfe/dinvn) #for fixed N and Te  

fit_p = np.dot(pressure_functions, coeff)

mid_hfe = 0.5*(hfe[1:] + hfe[:-1])
mid_cfe = 0.5*(cfe[1:] + cfe[:-1])

mid_h_gibbs = mid_hfe + mid_invn*hpressure
mid_c_gibbs = mid_cfe + mid_invn*cpressure
fit_c_gibbs = fit_cfe + invn*fit_p


# Find point of intersection for Figure 2
zoom_volume = 0.99
#plt.plot(fit_p, fit_c_gibbs - fit_p*zoom_volume, 'b:', label="fit crystal")

#Find pressure at point of intersection
def find_first_intersection(p1, g1, p2, g2): 
    for i in range(1,len(g1)-1):
        m1=(g1[i+1]-g1[i])/(p1[i+1]-p1[i])
        for j in range(1,len(g2)-1):
            m2=(g2[j+1]-g2[j])/(p2[j+1]-p2[j])
            if m1!=m2 :
                P_inter=(g2[j] - m2*p2[j] -g1[i] + m1*p1[i])/(m1-m2)
                if p1[i] < P_inter < p1[i+1] and p2[j] < P_inter < p2[j+1]:
                    g_inter=m1*P_inter+g1[i]-m1*p1[i]
                    if g1[i] < g_inter < g1[i+1] and g2[j] < g_inter < g2[j+1] :
                        return P_inter, g_inter


p_inter, g_inter = find_first_intersection(hpressure, mid_h_gibbs, cpressure, mid_c_gibbs)

pf_inter, gf_inter = find_first_intersection(hpressure, mid_h_gibbs, fit_p, fit_c_gibbs)


#Find homogeneous and crystal densities at p_inter
def find_densities(p_inter, pressure, invn):
    for i in range(1,len(pressure)-1): 
        if pressure[i] > p_inter :
            pressureabove=pressure[i]
            invnabove=invn[i]
            pressurebelow=pressure[i-1]
            invnbelow=invn[i-1]
            m=(pressureabove-pressurebelow)/(invnabove-invnbelow)
            invn_inter=invnabove-((pressureabove-p_inter)/m)
            return invn_inter
invnh=find_densities(p_inter, hpressure, mid_invn)
invnc=find_densities(p_inter, cpressure, mid_invn)


#-------------PLOTS---------------------------------------

# Plot Free Energy/atom vs 1/Reduced Density - Figure 1
#plt.plot(invn, fit_cfe, label="fit crystal free energy")
plt.plot(invn, hfe, 'red', label="Liquid")
plt.plot(invn, cfe, 'blue', label="Solid")
plt.xticks([invnc,invnh], [r'$\frac{1}{n_S}$', r'$\frac{1}{n_L}$'])
#plt.plot([invnc, invnc], [-163, -68.07], 'k-')
#plt.plot([invnh, invnh], [-163, -77.82], 'k-')
#plt.plot([invnc, invnh], [-68.09, -77.82], 'k-') 
#plt.plot([0.55, invnh], [-77.82 - ((invnh-0.55)*(-77.8-(-68.07))/(invnh-invnc)), -77.82], 'k-')
plt.title("Helmholtz Free Energy/atom vs 1/Reduced Density at Fixed kT=%g" % (kT))
plt.xlabel('Inverse Density')
#plt.xlabel(r'Inverse Reduced Density = $\frac{1}{n\sigma^3}$')
plt.ylabel('Helmholtz Free Energy/atom')
plt.legend()
plt.savefig('plot-pressure_Fig1.pdf')

plt.figure()

# Plot Gibbs Free Energy/atom vs Pressure with point of intersection   - Figure 2
#plt.plot(hpressure, mid_h_gibbs - hpressure*zoom_volume, 'red', label="Liquid")
#plt.plot(cpressure, mid_c_gibbs - cpressure*zoom_volume, 'blue', label="Solid")
#plt.xticks([p_inter], [r'$p_{Trans}$'])
#plt.yticks([g_inter- p_inter*zoom_volume], [r'$\mu_{Trans}$'])
#plt.plot([0, p_inter], [g_inter-p_inter*zoom_volume, g_inter-p_inter*zoom_volume], 'k-')
#plt.plot([p_inter, p_inter], [-1080, g_inter-p_inter*zoom_volume], 'k-')
#plt.plot(p_inter, g_inter - p_inter*zoom_volume, 'o', color='turquoise', markersize=10)
#plt.plot(pf_inter, gf_inter - pf_inter*zoom_volume, 'o', color='orange', markersize=10)
#plt.ylabel('Chemical Potential - pressure*(%g volume unit)' % zoom_volume)

plt.plot(hpressure, mid_h_gibbs, 'red', label="Liquid")
plt.plot(cpressure, mid_c_gibbs, 'blue', label="Solid")
plt.xticks([p_inter], [r'$p_{Trans}$'])
plt.yticks([g_inter], [r'$\mu_{Trans}$'])
plt.plot([0, p_inter], [g_inter, g_inter], 'k-')
plt.plot([p_inter, p_inter], [-137, g_inter], 'k-')
plt.plot(p_inter, g_inter, 'o', color='orange', markersize=10)
#plt.plot(pf_inter, gf_inter, 'o', color='orange', markersize=10)
plt.ylabel('Chemical Potential')

#plt.xlabel('Reduced Pressure')
plt.xlabel('Pressure')
plt.title("Chemical Potential vs Pressure at Fixed kT=%g" % (kT))
plt.legend()
plt.savefig('plot-pressure_Fig2.pdf')
#note: chemical potential = Gibbs Free Energy/atom

plt.figure()


# Plot Pressure vs 1/Reduced Density with point of intersection - Figure 3
plt.xticks([invnc,invnh], [r'$\frac{1}{n_S}$', r'$\frac{1}{n_L}$'])
plt.yticks([p_inter], [r'$p_{Trans}$'])
plt.plot([invnc, invnc], [0, p_inter], 'k-')
plt.plot([invnh, invnh], [0, p_inter], 'k-')
plt.plot([0, invnh], [p_inter, p_inter], 'k-')
plt.plot(invnh, p_inter, 'o', color='red', markersize=10) 
plt.plot(invnc, p_inter, 'o', color='blue', markersize=10) 
fit_p = np.dot(pressure_functions, coeff)
#plt.plot(invn, fit_p, 'g.-', label="fit crystal pressure")
plt.plot(mid_invn, hpressure,  'red', label="Liquid")
plt.plot(mid_invn, cpressure, 'blue', label="Solid")
plt.title("Reduced Pressure vs 1/Reduced Density at Fixed kT=%g" % (kT))
#plt.xlabel(r'Inverse Reduced Density = $\frac{1}{n\sigma^3}$')
#plt.ylabel(r'Reduced Pressure = $\frac{P\sigma^3}{\epsilon}$')
plt.xlabel('Inverse Density')
plt.ylabel('Pressure')
plt.legend()
plt.savefig('plot-pressure_Fig3.pdf')





