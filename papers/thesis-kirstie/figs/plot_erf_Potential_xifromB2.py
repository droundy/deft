#!/usr/bin/python3

#Run this program from the directory it is listed in
#with command ./plot_erf_Potential_xifromB2.py

from scipy import special
from scipy.special import iv
from scipy.special import erf
import scipy.integrate
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
import math

epsilon=1
sigma=1

twoR=2**(1.0/6)*sigma

PI=np.pi
Xi=[]
B2_WCA_list=[]
B2_erf_list=[]

def alpha(T) :
    return (2.0/(1+np.sqrt(T*np.log(2))))**(1.0/6)
    
def xi_Eric(T) :
	return alpha(T)/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/T)+np.log(2)))

#Compute B2_erf analytically (default method): 
def B2_erf_analytical(Xi,T) :   
    return np.pi/3*((alpha(T)**3 + 1.5*alpha(T)*Xi**2)*(1+erf(alpha(T)/Xi)) + 1/np.sqrt(np.pi)*(alpha(T)**2*Xi + Xi**3)*np.exp(-(alpha(T)/Xi)**2))
 
#Compute B2_wca numerically by evaluating the integral (default method):
def Vwca(r):
    return 4*((1/r)**12-(1/r)**6) +  1
    
def f_wca(r, T):
    return np.exp(-Vwca(r)/T) - 1

def B2_wca_numerical(T):
    rmax = 2.0**(1/6)
    r = np.linspace(0, rmax, 10000)
    f = f_wca(r, T)
    f[0] = 1
    return -0.5*(4*np.pi*r**2*r[1]*f).sum()  


def find_Xi(T):
    B2wca = B2_wca_numerical(T)   #use with B2 wca numerical method  (default)
    xi_lo = 0
    xi_hi = 1
    while xi_hi - xi_lo > 0.000001:
        xi_mid = 0.5*(xi_hi + xi_lo)
        if B2_erf_analytical(xi_mid, T) > B2wca:    #use with B2 erf analytical method (default)
            xi_hi = xi_mid
        else:
            xi_lo = xi_mid
    return xi_mid
    

T = np.linspace(0.2575, 10, 100)

for KbT in T:
    B2_WCA_at_T = B2_wca_numerical(KbT)     #use with B2 wca numerical method  (default)
    
    Xi_at_T = find_Xi(KbT)
                    
    Xi.append(Xi_at_T)
    B2_erf_list.append(B2_erf_analytical(Xi_at_T, KbT))    #use with B2 erf analytical method  (default)
    B2_WCA_list.append(B2_WCA_at_T)
    print(KbT, Xi_at_T)



#Plot Vwca and Verf - Figure 1
r=np.linspace(1, 1.12, 2000) #for erf plots #original SAVE
#r=np.linspace(.6, 1.1225, 2000)   #2R=1.1225  shows intersection 
#r=np.linspace(.95, 1.1225, 2000)   #2R=1.1225  good

i=0
line_colors = ['orange', 'g', 'b', 'r', 'm']
for KbT in [10] :
	
	#Plot WCA
    sigma_over_r_to_pow6 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)
    V_wca = 4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon
    plt.plot(r/sigma, V_wca/KbT, label='Vwca', color='black')
	
    #Plot erf with Xi from B2
    Verf=-KbT*np.log(.5*(special.erf((r-alpha(KbT))/find_Xi(KbT))+1))
    plt.plot(r/sigma, Verf/KbT, label='kT=%g' % (KbT), color=line_colors[3])

    #Plot erf with Eric's Xi 
    Verf=-KbT*np.log(.5*(special.erf((r-alpha(KbT))/xi_Eric(KbT))+1))
    plt.plot(r/sigma, Verf/KbT, linestyle='dashed', label='Eric kT=%g' % (KbT), color=line_colors[3])
       
    i = i + 1

plt.xlabel('r/$\sigma$')
plt.ylabel('V(r)/KT')
plt.title('Compare Vwca/kT and Verf/kT at kT=%g' % (KbT))
plt.legend()


plt.figure()



#Potential Plot at various temperatures - Figure 2
#r=np.linspace(.6, 1.1225, 2000)   #2R=1.1225  shows intersection
#r=np.linspace(.8, 1.1225, 2000)   #2R=1.1225  
r=np.linspace(.95, 1.1225, 2000)   #2R=1.1225   good
#r=np.linspace(1, 1.12, 2000)   #2R=1.1225   matches old figure

#Plot WCA
sigma_over_r_to_pow6 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)
V_wca = 4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon
plt.plot(r/sigma, V_wca/epsilon, label='Vwca', color='black')

i=0
line_colors = ['orange', 'g', 'b', 'r', 'm']
for KbT in [0.5, 1, 2, 10, 40] :
	
    #Plot erf with Xi from B2
    Verf=-KbT*np.log(.5*(special.erf((r-alpha(KbT))/find_Xi(KbT))+1))
    plt.plot(r/sigma, Verf/epsilon, label='kT=%g' % (KbT), color=line_colors[i])

    #Plot erf with Eric's Xi 
    Verf=-KbT*np.log(.5*(special.erf((r-alpha(KbT))/xi_Eric(KbT))+1))
    #plt.plot(r/sigma, Verf/epsilon, linestyle='dashed', label='Eric kT=%g' % (KbT), color=line_colors[i])
       
    i = i + 1
    
plt.xlabel('r/$\sigma$')
plt.ylabel('V(r)/$\epsilon$')
plt.title('Compare Vwca and Verf')
plt.legend()


plt.show()


