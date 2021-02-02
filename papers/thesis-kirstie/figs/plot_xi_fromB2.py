#!/usr/bin/python3

#Run this program from the directory it is listed in
#with command ./plot_xi_fromB2.py

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
	return alpha(T)/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/T)+np.log(2)))    #old xi found from matching slope (Eric's method)

#Compute B2_erf analytically (default method): 
def B2_erf_analytical(Xi,T) :   
    return np.pi/3*((alpha(T)**3 + 1.5*alpha(T)*Xi**2)*(1+erf(alpha(T)/Xi)) + 1/np.sqrt(np.pi)*(alpha(T)**2*Xi + Xi**3)*np.exp(-(alpha(T)/Xi)**2))
    
#NEW Compute B2_erf analytically  (experimental): YUK!
def new_B2_erf_analytical(Xi,T) :
    new_alpha=pow(2,1.0/6)-2*Xi
#    new_alpha=pow(2,1.0/6)-(1.5+0.3*np.log(T))*Xi
    return np.pi/3*((new_alpha**3 + 1.5*new_alpha*Xi**2)*(1+erf(new_alpha/Xi)) + 1/np.sqrt(np.pi)*(new_alpha**2*Xi + Xi**3)*np.exp(-(new_alpha/Xi)**2))    
    
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
        #if new_B2_erf_analytical(xi_mid, T) > B2wca:    #use with B2 erf analytical method (experimental) YUK!
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



#Plot B2_WCA vs T - Figure 1
plt.plot(T, B2_WCA_list, label = 'B2_WCA')
plt.plot(T, B2_erf_list, label = 'B2_erf')
plt.xlabel('KbT')
plt.ylabel('B2')
plt.xlim(0,T.max())
plt.ylim(bottom=0)
plt.title('B2 vs Temp')
#plt.legend()
plt.savefig("B2_vs_T_with_xi_fromB2.pdf")

plt.figure()



#Plot Xi vs T - Figure 2
plt.plot(T, Xi, label = 'Xi_B2', color='blue')
#plt.plot(T, xi_Eric(T), label='Krebs', color='orange')
#plt.plot(T, (0.12**2*(T-0.2575))**(1/2.1), label='crazy fit')
plt.xlabel('KbT')
plt.ylabel('Xi')
plt.xlim(0,T.max())
plt.ylim(bottom=0)
plt.title('Xi vs Temp')
#plt.legend()
plt.savefig("xi_fromB2.pdf")


plt.figure()



#Plot df/dr and w2*w2 for comparison - Figure 3
r=np.linspace(.6, 1.1225, 2000) #Use this for full range of temperatures
#r=np.linspace(0.9, 1.20, 2000) #Use to match Eric's diagram
#r=np.linspace(0.7, 1.25, 2000) #

KbT=2 #this is actually kBT/epsilon with epsilon=1

print(alpha(KbT))

w2_w2_Eric=(1.0/(xi_Eric(KbT)*np.sqrt(PI)))*np.exp(-pow((r-alpha(KbT))/xi_Eric(KbT),2))  #convolution of w2 with itself
w2_w2_B2=(1.0/(find_Xi(KbT)*np.sqrt(PI)))*np.exp(-pow((r-alpha(KbT))/find_Xi(KbT),2))  #convolution of w2 with itself  (use this with Xi from B2)

df_dr= np.exp(-(1/KbT)*(4*epsilon*(pow(sigma,12)*pow(r,-12)-pow(sigma,6)*pow(r,-6))+epsilon))*((1/KbT)*4*epsilon*(12*pow(sigma,12)*pow(r,-13)-6*pow(sigma,6)*pow(r,-7))) #derivative of the mayer function from the WCA potential

plt.plot(r/alpha(KbT), df_dr, label='WCA df_dr', color='black')
#plt.plot(r/sigma, w2_w2_Eric, label='Eric w2_w2', linestyle='dashed', color='orange')
plt.plot(r/alpha(KbT), w2_w2_B2, label='w2*w2', color='blue') #with Xi from matching B2

plt.xlim(left=0.6,right=1.12)
plt.ylim(bottom=0,top=8.2)

#plt.xticks([0], [' '])
plt.xticks([1], [])
plt.yticks([0], [' '])


#plt.xlabel('r/$\sigma$')
#plt.xlabel('r/alpha')
#plt.ylabel('df(r)/dr, w2*w2')
#plt.title('Compare match with Eq12 for low density kT=%g' % (KbT))
plt.legend()
plt.savefig("df_dr_with_xi_fromB2.pdf")


#plt.show()


