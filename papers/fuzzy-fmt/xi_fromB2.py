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
#T=[]
Xi=[]
B2_WCA_list=[]
B2_erf_list=[]

def alpha(KbT) :
    return pow(np.sqrt((2.0/(1+np.sqrt(KbT*np.log(2))))), 1.0/3)

def B2_WCA_integrand(r, KbT) :
    V_WCA=4*(pow(r,-12) - pow(r,-6)) +1
    WCA_mayer_function=np.exp(-V_WCA/KbT)-1
    return (-1.0/2)*4*np.pi*WCA_mayer_function*r*r
    
def B2_erf_integrand(r, KbT, pXi) :
    erf_mayer_function=(1.0/2)*(erf((r-alpha(KbT))/(pXi/np.sqrt(2)))-1)
    return (-1.0/2)*(4*np.pi)*(erf_mayer_function)*r*r

def B2_wca(T):
    upper_limit=sigma*pow(2,1.0/6)
    return quad(B2_WCA_integrand, 0, upper_limit, args=(KbT))  #integrate from 0 to less than 2R

def B2_erf(T, Xi):
    return quad(B2_erf_integrand, 0, np.inf, args=(T, Xi))

def find_Xi(T):
    B2wca = B2_wca(KbT)[0]
    xi_lo = 0
    xi_hi = 1
    while xi_hi - xi_lo > 0.000001:
        xi_mid = 0.5*(xi_hi + xi_lo)
        if B2_erf(T, xi_mid)[0] > B2wca:
            xi_hi = xi_mid
        else:
            xi_lo = xi_mid
    return xi_mid

T = np.linspace(0.2575, 10, 100)
for KbT in T:

    B2_WCA_at_T = B2_wca(KbT)[0]
    
    Xi_at_T = find_Xi(KbT)
                    
    #print "B2_diff=", B2_diff, "Xi_at_T=", Xi_at_T, "KbT=", KbT
    Xi.append(Xi_at_T)
    B2_erf_list.append(B2_erf(KbT, Xi_at_T)[0])
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

alpha=np.cbrt(np.sqrt((2/(1+np.sqrt(T*np.log(2))))))
Xi_old = alpha/(6*np.sqrt(np.pi)*(np.sqrt(np.log(2)/T) + np.log(2)))

##Plot xi
plt.plot(T, Xi, label = 'Xi')
plt.plot(T, Xi_old, label='Krebs')
#plt.plot(T, (0.09**2*(T-0.2575))**(1/2), label='crazy fit')
plt.xlabel('KbT')
plt.ylabel('Xi')
plt.title('Xi vs Temp')
plt.legend()

plt.show()


# #-----CHECK--------------------
# T_ch=[]
# B2_WCA_ch=[]
# B2_erf_ch=[]

# pXi=0.171
# for KbT in np.arange(0.05, 300, 0.15):
    # B2_erf_at_T=quad(B2_erf_integrand, 0, np.inf, args=(KbT, pXi))
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
# plt.title('B2 vs T at Xi=%g' % (pXi))
# plt.legend()

# plt.show()
