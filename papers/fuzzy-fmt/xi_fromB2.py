#!/usr/bin/python3

#Run this program from the directory it is listed in
#with command ./xi_fromB2.py

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

PI=np.pi
Xi=[]
B2_WCA_list=[]
B2_erf_list=[]

def alpha(T) :
    return (2.0/(1+np.sqrt(T*np.log(2))))**(1.0/6)

#Compute B2_erf analytically (default method): 
def B2_erf_analytical(Xi,T) :   
    return np.pi/3*((alpha(T)**3 + 1.5*alpha(T)*Xi**2)*(1+erf(alpha(T)/Xi)) + 1/np.sqrt(np.pi)*(alpha(T)**2*Xi + Xi**3)*np.exp(-(alpha(T)/Xi)**2))
    
#Compute B2_erf numerically by evaluating the integral:    
def f_erf(r, Xi, T):
    return (0.5)*(erf((r-alpha(T))/Xi)-1)

def B2_erf_numerical(Xi,T):
    r = np.linspace(0, 10*alpha(T), 10000)
    f = f_erf(r, Xi, T)
    f[0] = 1
    return -0.5*(4*np.pi*r**2*r[1]*f).sum() 
    
#Compute B2_erf by using python quad function:
def B2_erf_quad_integrand(r, Xi, T) :
    erf_mayer_function=(1.0/2)*(erf((r-alpha(T))/Xi)-1)
    return (-1.0/2)*(4*np.pi)*(erf_mayer_function)*r*r
    
def B2_erf_quad(Xi, T):
    return quad(B2_erf_quad_integrand, 0, np.inf, args=(Xi, T))


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
        
#Compute V2_wca by using python quad function:
def B2_WCA_quad_integrand(r, T) :
    V_WCA=4*(pow(r,-12) - pow(r,-6)) +1
    WCA_mayer_function=np.exp(-V_WCA/T)-1
    return (-1.0/2)*4*np.pi*WCA_mayer_function*r*r

def B2_wca_quad(T):
    rmax=sigma*pow(2,1.0/6)
    return quad(B2_WCA_quad_integrand, 0, rmax, args=(T))  #integrate from 0 to less than 2R


def find_Xi(T):
    #B2wca = B2_wca_quad(T)[0]   #use with B2 wca quad method
    B2wca = B2_wca_numerical(T)   #use with B2 wca numerical method  (default)
    xi_lo = 0
    xi_hi = 1
    while xi_hi - xi_lo > 0.000001:
        xi_mid = 0.5*(xi_hi + xi_lo)
        #if B2_erf_quad(xi_mid, T)[0] > B2wca:      #use with B2 erf quad method
        if B2_erf_analytical(xi_mid, T) > B2wca:    #use with B2 erf analytical method (default)
        #if B2_erf_numerical(xi_mid, T) > B2wca:    #use with B2 erf numerical method
            xi_hi = xi_mid
        else:
            xi_lo = xi_mid
    return xi_mid
    

T = np.linspace(0.2575, 10, 100)

for KbT in T:
    #B2_WCA_at_T = B2_wca_quad(KbT)[0]      #use with B2 wca quad method
    B2_WCA_at_T = B2_wca_numerical(KbT)     #use with B2 wca numerical method  (default)
    
    Xi_at_T = find_Xi(KbT)
                    
    Xi.append(Xi_at_T)
    #B2_erf_list.append(B2_erf_quad(Xi_at_T, KbT)[0])      #use with B2 erf quad method
    B2_erf_list.append(B2_erf_analytical(Xi_at_T, KbT))    #use with B2 erf analytical method  (default)
    #B2_erf_list.append(B2_erf_numerical(Xi_at_T, KbT))    #use with B2 erf numerical method
    B2_WCA_list.append(B2_WCA_at_T)
    print(KbT, Xi_at_T)
   
#Plot B2_WCA vs T
plt.plot(T, B2_WCA_list, label = 'B2_WCA')
plt.plot(T, B2_erf_list, label = 'B2_erf')
plt.xlabel('KbT')
plt.ylabel('B2')
plt.title('B2 vs Temp')
plt.legend()

plt.figure()


Xi_old = alpha(T)/(6*np.sqrt(np.pi)*(np.sqrt(np.log(2)/T) + np.log(2)))

##Plot xi
plt.plot(T, Xi, label = 'Xi')
plt.plot(T, Xi_old, label='Krebs')
plt.plot(T, (0.12**2*(T-0.2575))**(1/2.1), label='crazy fit')
plt.xlabel('KbT')
plt.ylabel('Xi')
plt.title('Xi vs Temp')
plt.legend()

plt.show()


# #-----CHECK--------------------
# T_ch=[]
# B2_WCA_ch=[]
# B2_erf_ch=[]

# Xi=0.171
# for KbT in np.arange(0.05, 300, 0.15):
    # B2_erf_at_T=quad(B2_erf_integrand, 0, np.inf, args=(Xi, KbT))
    # B2_WCA_at_T=quad(B2_WCA_integrand, 0, 2/1.781797435, args=(KbT))
    # B2_erf_ch.append(B2_erf_at_T[0])
    # B2_WCA_ch.append(B2_WCA_at_T[0])
    # T_ch.append(KbT)
    
# #plot B2 vs T
# # plt.plot(T_ch, B2_WCA_ch, label = 'B2_WCA_ch')
# # plt.plot(T_ch, B2_erf_ch, label = 'B2_erf_ch')
# plt.plot(T_ch, np.array(B2_erf_ch)-np.array(B2_WCA_ch), label = 'B2_WCA_ch')
# plt.axhline(0)
# plt.xlabel('KbT')
# plt.ylabel('B2')
# plt.title('B2 vs T at Xi=%g' % (Xi))
# plt.legend()

# plt.show()

