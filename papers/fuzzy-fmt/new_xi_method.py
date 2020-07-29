#!/usr/bin/python3

#Run this program from the directory it is listed in
#with command ./new_xi_method.py

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
    print(KbT, Xi_at_T)


#Range of "magic #" for good results is 1.4 to 1.6, my choice 1.5 or 1.43
Xi_old = alpha(T)/(6*np.sqrt(np.pi)*(np.sqrt(np.log(2)/T) + np.log(2)))
Xi_new_2 = (pow(2,1.0/6)-alpha(T))/2.0 #not good
Xi_new_1_6 = (pow(2,1.0/6)-alpha(T))/1.6  #SAVE 1.6 works ok
Xi_new_1_5 = (pow(2,1.0/6)-alpha(T))/1.5  #save 1.5 better - crosses Xi_B2 twice!
#Xi_new_1_45 = (pow(2,1.0/6)-alpha(T))/1.45  #also good!
Xi_new_1_44 = (pow(2,1.0/6)-alpha(T))/1.44  #also good!
Xi_new_1_43 = (pow(2,1.0/6)-alpha(T))/1.43 #close fit for 1.1<T<5
#Xi_new_1_42 = (pow(2,1.0/6)-alpha(T))/1.42 #close fit for 1.3<T<4 staring to go over Xi_B2
#Xi_new_1_41 = (pow(2,1.0/6)-alpha(T))/1.41 #close fit for 1.4<T<3.5 staring to go over Xi_B2
Xi_new_1_4 = (pow(2,1.0/6)-alpha(T))/1.40 #SAVE close fit for 1.6<T<2.9 staring to go over Xi_B2
Xi_new_ln=(pow(2,1.0/6)-alpha(T))/(1.5+0.3*np.log(T)) #Xi_new_ln
Xi_new_ln10=(pow(2,1.0/6)-alpha(T))/(1.5+0.3*np.log(10)) #Xi_new_ln10   #took pictures of this!
#Xi_new_ln=(pow(2,1.0/6)-alpha(T))/(1.5+0.3) #Xi_new
print(T)
#Xi_df_dr=(T/(4*np.sqrt(PI)))*(1.0/(12*pow(alpha(T),-13)-6*pow(alpha(T),-7)))*np.exp((1/T)*(4*pow(alpha(T),-12)-4*pow(alpha(T),-6)+1))  #Xi_df_dr
#OTHER: Xi_new = pow(2,-5.0/6)-(1.0/2.15)*alpha(T)  #not correct equation, but interesting


##Plot xi
plt.plot(T, Xi, label = 'Xi_B2', color='blue')
plt.plot(T, Xi_old, label='Krebs', color='orange')
#plt.plot(T, Xi_new, label='Xi_2pt new')
#plt.plot(T, Xi_new_2, label='Xi_2pt 2', color='green' )
#plt.plot(T, Xi_new_1_6, label='Xi_2pt 1.6')
#plt.plot(T, Xi_new_1_5, label='Xi_2pt 1.5')
#plt.plot(T, Xi_new_1_43, label='Xi_2pt 1.43', color='yellow'  )
#plt.plot(T, Xi_new_1_4, label='Xi_2pt 1.4')
plt.plot(T, Xi_new_ln, label='Xi ln', color='red' )
plt.plot(T, Xi_new_ln10, label='Xi ln10', color='purple' )
plt.plot(T, Xi_new_1_43, label='Xi_2pt', color='yellow' )
#plt.plot(T, Xi_df_dr, label='Xi_df_dr', color='purple' )
#plt.plot(T, (0.12**2*(T-0.2575))**(1/2.1), label='crazy fit')
plt.xlabel('KbT')
plt.ylabel('Xi')
plt.title('Xi vs Temp')
plt.legend()

plt.figure()

#Plots---------------
KbT=1 #this is actually kBT/epsilon with epsilon=1


zeta_B2 = find_Xi(KbT)  #Xi from B2
zeta_Eric = alpha(KbT)/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/KbT)+np.log(2)))
zeta_ln=(pow(2,1.0/6)-alpha(KbT))/(1.5+0.3*np.log(KbT)) #Xi_new_ln
zeta_ln10=(pow(2,1.0/6)-alpha(KbT))/(1.5+0.3*np.log(10))
zeta_2pt=(pow(2,1.0/6)-alpha(KbT))/1.43 #Xi_new_2pt
zeta_df_dr=(KbT/(4*np.sqrt(PI)))*pow(12*pow(alpha(KbT),-13)-6*pow(alpha(KbT),-7),-1.0)*np.exp((1/KbT)*(4*pow(alpha(KbT),-12)-4*pow(alpha(KbT),-6)+1))  #Xi_df_dr

alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))

#Potential Plots
#r=np.linspace(1, 1.12, 2000) #for erf plots #original SAVE
r=np.linspace(.5, 1.1225, 2000)   #2R=1.1225
#r=np.linspace(.5, 3, 2000)

#Plot WCA
sigma_over_r_to_pow6 = (sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)
V = 4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon
#plt.plot(r/sigma,V/epsilon, label='Vwca', color='black')
plt.plot(r/sigma,V/KbT, label='Vwca V/T', color='black')

#Plot erf with Xi from B2  
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_B2)+1))
#plt.plot(r/sigma,Verf/epsilon, label='B2', color='blue')
plt.plot(r/sigma,Verf/KbT, label='B2 V/T', color='blue')

#Plot erf with Xi_ln 
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_ln)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi ln', color='red')
plt.plot(r/sigma,Verf/KbT, label='Xi ln V/T', color='red')

#Plot erf with Xi_ln10 
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_ln10)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi ln', color='red')
plt.plot(r/sigma,Verf/KbT, label='Xi ln10 V/T', color='purple')

#Plot erf with Eric's Xi 
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_Eric)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Eric', color='orange')
plt.plot(r/sigma,Verf/KbT, label='Eric V/T', color='orange')

#Plot erf with Xi_2pt 1.43
#zeta=(pow(2,1.0/6)-alpha)/1.43 #Xi_new_1_43
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_2pt)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi 1.43', color='yellow')
plt.plot(r/sigma,Verf/KbT, label='Xi 1.43 V/T', color='yellow')

#Plot erf with Xi_2pt 2
zeta=(pow(2,1.0/6)-alpha)/2 #Xi_new_2
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi 2', color='green')

#Plot erf with Xi_2pt 3 
zeta=(pow(2,1.0/6)-alpha)/3 #Xi_new_3
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi 3', color='yellow')

#Plot erf with Xi_df_dr 
Verf=-KbT*np.log(.5*(special.erf((r-alpha)/zeta_df_dr)+1))
#plt.plot(r/sigma,Verf/epsilon, label='Xi ln', color='purple')

plt.xlabel('r/$\sigma$')
#plt.ylabel('V(r)/$\epsilon$')
plt.ylabel('V(r)/KT')
plt.title('Error Function Potential Verf kT=%g' % (KbT))
plt.legend()

plt.figure()

#Plot df/dr and w2*w2 for comparison---------------------
r=np.linspace(.3, 1.1225, 2000) #Use this for full range of temperatures
#r=np.linspace(0.9, 1.20, 2000) #Use to match Eric's diagram

w2_w2_Eric=(1.0/(zeta_Eric*np.sqrt(PI)))*np.exp(-pow((r-alpha)/zeta_Eric,2))  #convolution of w2 with itself
w2_w2_B2=(1.0/(zeta_B2*np.sqrt(PI)))*np.exp(-pow((r-alpha)/zeta_B2,2))  #convolution of w2 with itself  (use this with Xi from B2)
w2_w2_ln=(1.0/(zeta_ln*np.sqrt(PI)))*np.exp(-pow((r-alpha)/zeta_ln,2))
w2_w2_ln10=(1.0/(zeta_ln10*np.sqrt(PI)))*np.exp(-pow((r-alpha)/zeta_ln10,2))
w2_w2_2pt=(1.0/(zeta_2pt*np.sqrt(PI)))*np.exp(-pow((r-alpha)/zeta_2pt,2))

df_dr= np.exp(-(1/KbT)*(4*epsilon*(pow(sigma,12)*pow(r,-12)-pow(sigma,6)*pow(r,-6))+epsilon))*((1/KbT)*4*epsilon*(12*pow(sigma,12)*pow(r,-13)-6*pow(sigma,6)*pow(r,-7))) #derivative of the mayer function from the WCA potential

plt.plot(r/sigma, df_dr, label='WCA df_dr')
plt.plot(r/sigma, w2_w2_Eric, label='Eric w2_w2', linestyle='dashed', color='orange')
plt.plot(r/sigma, w2_w2_B2, label='B2 w2_w2', linestyle='dashed', color='blue')
plt.plot(r/sigma, w2_w2_ln, label='ln w2_w2', linestyle='dashed', color='red')
plt.plot(r/sigma, w2_w2_ln10, label='ln10 w2_w2', linestyle='dashed', color='purple')
plt.plot(r/sigma, w2_w2_2pt, label='2pt w2_w2', linestyle='dashed', color='yellow')
plt.xlabel('r/$\sigma$')
plt.ylabel('df(r)/dr, w2*w2')
plt.title('Compare match with Eq12 for low density kT=%g' % (KbT))
plt.legend()


plt.show()

print('alpha=',alpha)




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
