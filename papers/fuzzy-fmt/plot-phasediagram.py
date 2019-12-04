#!/usr/bin/python2

#This program produces temperature vs density, and pressure vs temperature phase diagrams
#from data stored in *best.dat (or *best_tensor.dat) data files generated by figs/new-melting.cpp
#and found in deft/papers/fuzzy-fmt/data/phase-diagram  (edit later - currently files in newdata/phase-diagram and newdata_tensor/phasediagram)

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./plot-phasediagram.py [Optional: --tensor]  

from __future__ import print_function, division

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse
import sys

parser = argparse.ArgumentParser("Plots phase diagrams p vs T and T vs n. Plots p-vs-T, p-vs-n, and T-vs-n.")
parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')

args=parser.parse_args()

p_at_freezing = []  #pressure at freezing (intersection point between homogeneous and crystal plots)
n_homogeneous_at_freezing =[]
n_crystal_at_freezing = []
kT_homogeneous_at_freezing = []
kT_crystal_at_freezing = []
kT_in_plot = []

kT_data = []
density_data = []   #index corresponds to kT
pressure_data = []  #index corresponds to kT

#for kT in np.arange(0.1, 1.15, 0.05):   #data files with these temperatures will be plotted
#for kT in np.arange(0.1, 2.05, 0.05):  #original
#for kT in np.arange(0.4, 2.05, 0.05):   # new normal
for kT in (1, 2, 4, 6, 8, 10, 12, 14, 16, 18):
#for kT in np.arange(0.1, 1.05, 0.05):   #data files with these temperatures will be plotted  DEBUG
										#values above and below this range do not currrently work   DEBUG
   
   n = []
   invn = []
   hfe = []
   cfe = []

   if args.tensor :
     #files = sorted(list(glob.glob('data/phase-diagram/kT%.3f_n*_best_tensor.dat' % kT)))
     files = sorted(list(glob.glob('newdata_tensor/phase-diagram3/kT%.3f_n*_best_tensor.dat' % kT)))    #remove 2 at the end of phase-diagram when done comparing new data
     
   else :
      files = sorted(list(glob.glob('newdata/phase-diagram3/kT%.3f_n*_best.dat' % kT)))    #remove 2 at the end of phase-diagram when done comparing new data
      #files = sorted(list(glob.glob('data/phase-diagram/kT%.3f_n*_best.dat' % kT)))
      #files = sorted(list(glob.glob('crystallization/kT%.3f_n*_best.dat' % kT)))

   if len(files) == 0:
	   continue
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


   #Find pressure at point of intersection
   def find_first_intersection(p1, g1, p2, g2): 
      for i in range(1,len(g1)-1):
         m1=(g1[i+1]-g1[i])/(p1[i+1]-p1[i])
         for j in range(1,len(g2)-1):
               m2=(g2[j+1]-g2[j])/(p2[j+1]-p2[j])
               #print(m2) #debug ASK!
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
  
   p_at_freezing.append(p_inter)   
   n_homogeneous_at_freezing.append(1/invnh)
   n_crystal_at_freezing.append(1/invnc)
   

   # compute the actual physical pressure as a function of density, and skip over coexistence
   actual_pressure = []
   actual_density = []
   for i in range(len(mid_invn)):
      if hpressure[i] >= p_inter:
         break # if the pressure is too high, then we should just stop, since we have left the fluid
      actual_pressure.append(hpressure[i])
      actual_density.append(1/mid_invn[i])
   actual_pressure.append(p_inter)
   actual_density.append(1/invnh)
   actual_pressure.append(p_inter)
   actual_density.append(1/invnc)
   for i in range(len(mid_invn)):
      if cpressure[i] < 0 and mid_invn[i] <= invnc:
         break # when the pressure is negative, we know we are in the crazy part where our dft fails.
      if cpressure[i] > p_inter:
         actual_pressure.append(cpressure[i])
         actual_density.append(1/mid_invn[i])
   actual_pressure = np.array(actual_pressure)
   actual_density = np.array(actual_density)

   #print (kT, p_inter, 1/invnh, 1/invnc)   #Use >> phase_diagram_data.dat (or phase_diagram_data-tensor.dat) to store data for reference
   
   kT_data.append(kT)  #holds all values of kT in a list
   density_data.append(actual_density)
   pressure_data.append(actual_pressure)

n_homogeneous_at_freezing = np.array(n_homogeneous_at_freezing)
n_crystal_at_freezing = np.array(n_crystal_at_freezing)
p_at_freezing = np.array(p_at_freezing)

plt.figure('T-vs-n at fixed P')
plt.fill_betweenx(kT_data, n_homogeneous_at_freezing, n_crystal_at_freezing, color='#eeeeee') 

#Plot T vs n  at constant P 
for p in [2,5,10,20]:
   n_mid_at_p_list = []
   kT_at_p_list = []
   for i in range(0, len(kT_data)) :  #number of temperatures kT
      for j in range(0, len(density_data[i])-1) :  #number of elements of n at some kT
         if pressure_data[i][j] < p < pressure_data[i][j+1] :
            phi = pressure_data[i][j+1]
            plo = pressure_data[i][j]
            nhi = density_data[i][j+1]
            nlo = density_data[i][j]
            n_mid_at_p_list.append((nlo*(phi - p) + nhi*(p - plo))/(phi - plo))
            kT_at_p_list.append(kT_data[i])

   plt.plot(n_mid_at_p_list, kT_at_p_list, '.-', label= 'P=%g' % p)
plt.title("Temperature vs Number Density at fixed Pressure")
plt.legend(loc='best')
plt.xlabel('Number Density')
plt.ylabel('Temperature')

# - OR - uncomment the plot you want
   
   #Plot n vs T  at constant P
   #plt.plot(kT_at_p_list, n_mid_at_p_list, '.-', label= 'P=%g' % p)
#plt.title("Number Density vs Temperature at fixed Pressure")
#plt.legend(loc='best')
#plt.ylabel('Number Density')
#plt.xlabel('Temperature') 
   
plt.figure('p-vs-n at fixed T')

plt.fill_betweenx(p_at_freezing, n_homogeneous_at_freezing, n_crystal_at_freezing, color='#eeeeee') 
for i in range(len(kT_data)):
   if kT_data[i] in [0.1, 0.2, 0.5, 1.0] or True:
      #Plot P vs n  at constant kT
      plt.plot(density_data[i], pressure_data[i], label= 'kT=%g' % kT_data[i])
plt.title("Pressure vs Number Density at kT")
plt.legend(loc='best')
#plt.ylim(0, 26)
plt.ylim(0, 500)
#plt.xlim(0, 1.1)
plt.xlim(0, 1.8)
plt.xlabel('Number Density')
plt.ylabel('Pressure')

plt.figure('p-vs-V at fixed T')

#Plot P vs 1/n (or V) at constant kT
plt.fill_betweenx(p_at_freezing, 1/n_homogeneous_at_freezing, 1/n_crystal_at_freezing, color='#eeeeee') 
for i in range(len(kT_data)):
   if kT_data[i] in [0.1, 0.2, 0.5, 1.0] or True:
      plt.plot(1/density_data[i], pressure_data[i], label= 'kT=%g' % kT_data[i])
plt.title("Pressure vs volume at kT")
plt.legend(loc='best')
plt.ylim(0, 26)
plt.xlim(0.95, 1.6)
plt.xlabel('Volume per atom')
plt.ylabel('Pressure')

plt.figure('p-vs-T at fixed n')

#--------------NEW
#Plot P vs T  at constant n 
#for n in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]:    #densities to show on the plot
for n in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1]:  #densities to show on the plot
   p_mid_at_n_list = []
   kT_at_n_list = []
   for i in range(0, len(kT_data)) :  #number of temperatures kT
      for j in range(0, len(pressure_data[i])-1) :  #number of elements of P at some kT
         if density_data[i][j] < n < density_data[i][j+1] :
            phi = pressure_data[i][j+1]
            plo = pressure_data[i][j]
            nhi = density_data[i][j+1]
            nlo = density_data[i][j]
            p_mid_at_n_list.append((plo*(nhi - n) + phi*(n - nlo))/(nhi - nlo))
            kT_at_n_list.append(kT_data[i])

   plt.plot(kT_at_n_list, p_mid_at_n_list, '.-', label= 'n=%g' % n)
plt.title("Pressure vs Temperature at fixed n")
plt.legend(loc='best')
plt.ylabel('Pressure')
plt.xlabel('Temperature')

# - OR - uncomment the plot you want
   ##Plot T vs P  at constant n 
   #plt.plot(kT_at_n_list, p_mid_at_n_list, '.-', label= 'n=%g' % n)
#plt.title("Temperature vs Pressure at fixed n")
#plt.legend(loc='best')
#plt.xlabel('Pressure')
#plt.ylabel('Temperature')
#--------------end NEW

plt.figure('Phase Diagram of T vs n')

#Temperature vs Density Phase Diagram
plt.plot(n_homogeneous_at_freezing, kT_data, label='liquid', color='red')
plt.plot(n_crystal_at_freezing, kT_data, label='solid', color='blue')
plt.fill_betweenx(kT_data, .1, n_homogeneous_at_freezing, color='red')       
plt.fill_betweenx(kT_data, n_homogeneous_at_freezing, n_crystal_at_freezing, color='gray') 
#plt.fill_betweenx(kT_data, n_crystal_at_freezing, 1.6, color='blue')
plt.fill_betweenx(kT_data, n_crystal_at_freezing, 1.8, color='blue')           
plt.title("Temperature vs Number Density")
#plt.legend(loc='best')
plt.xlabel('Number Density')
plt.ylabel('Temperature')

#plt.plot([0.88, 0.90, 0.91, 0.92, 1.04, 1.12],[0.7, 0.8, 0.9, 1.0, 2.0, 3.0], label='chris_l', color='green')
plt.plot([0.88, 0.90, 0.91, 0.92, 1.04, 1.12, 1.24, 1.44],[0.7, 0.8, 0.9, 1.0, 2.0, 3,5,10], label='chris_l', color='green')
#plt.plot([0.96, 0.98, 0.99, 1.00, 1.11, 1.19],[0.7, 0.8, 0.9, 1.0, 2.0, 3.0], label='chris_s', color='green')
plt.plot([0.96, 0.98, 0.99, 1.00, 1.11, 1.19, 1.31, 1.51],[0.7, 0.8, 0.9, 1.0, 2.0, 3, 5, 10], label='chris_s', color='green')
plt.legend()

plt.figure('Phase Diagram of P vs T')

##Pressure vs Temperature Phase Diagram
plt.fill_between(kT_data, 0, p_at_freezing, color='red')      
#plt.fill_between(kT_data, p_at_freezing, 50, color='blue')    #FIX - change 30
plt.fill_between(kT_data, p_at_freezing, 1500, color='blue')    
plt.plot(kT_data, p_at_freezing, color='black')
#plt.ylim(0, 40)
#plt.xlim(kT_data.min(), kT_data.max())     #FIX!  
plt.title("Pressure vs Temperature")
plt.xlabel('Temperature')
plt.ylabel('Pressure')

#plt.plot([0.7, 0.8,0.9,1.0,2.0,3.0], [6.24, 7.62, 8.78, 9.99, 25.5,43.8], label='chris_l', color='green')
plt.plot([0.7, 0.8,0.9,1.0,2.0, 3, 5, 10], [6.24, 7.62, 8.78, 9.99, 25.5,43.8, 85.6, 210], label='chris_l', color='green')
plt.legend()

plt.show()

