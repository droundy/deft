#!/usr/bin/python2

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
T=[]
Xi=[]
B2_WCA=[]
B2_erf=[]


def B2_WCA_integrand(r, KbT) :
    V_WCA=4*(pow(r,-12) - pow(r,-6)) +1
    WCA_mayer_function=np.exp(-V_WCA/KbT)-1
    return (-1.0/2)*4*np.pi*WCA_mayer_function*r*r
    
def B2_erf_integrand(r, KbT, pXi) :
    alpha=np.cbrt(np.sqrt((2/(1+np.sqrt(KbT*np.log(2))))))
    erf_mayer_function=(1.0/2)*(erf((r-alpha)/pXi)-1)
    return (-1.0/2)*(4*np.pi)*(erf_mayer_function)*r*r
    
upper_limit=2.0*sigma*pow(2,-5.0/6)   #integrate from 0 to less than 2R  where R=sigma*2^(-5/6)   FIX?
print upper_limit


for KbT in np.arange(0.5, 200, 0.5):

    B2_WCA_int = quad(B2_WCA_integrand, 0, upper_limit, args=(KbT))  #integrate from 0 to less than 2R
    
    alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))
    
    #METHOD 1: -----------------------------
    # #B2_diff=10 
    # for pXi in np.arange(0.01, 0.2, 0.01) :   #pXi stands for "proposed Xi"
    # #The real Xi at a particular KbT is found when B2_WCA=B2_erf at that particular KbT
        # pB2_erf_int=quad(B2_erf_integrand, 0, np.inf, args=(KbT, pXi))
        # pB2_diff=pB2_erf_int[0]-B2_WCA_int[0]
        # print(pB2_diff)
        # #if math.copysign(1, pB2_diff) =
        # if abs(pB2_diff) < B2_diff :
            # B2_diff=pB2_diff
            # B2_WCA_at_T=B2_WCA_int[0]
            # B2_erf_at_T=pB2_erf_int[0]
            # Xi_at_T=pXi 
            #-------------------------------   
                       
    #METHOD 2: -----------------------------
    pB2_erf_int=quad(B2_erf_integrand, 0, np.inf, args=(KbT, 0.01))
    pB2_diff=pB2_erf_int[0]-B2_WCA_int[0]
    if pB2_diff > 0 :
        sign=1
    elif pB2_diff < 0 :
        sign=0
    for pXi in np.arange(0.011, 0.2, 0.001) :   #pXi stands for "proposed Xi"
    #The real Xi at a particular KbT is found when B2_WCA=B2_erf at that particular KbT
        pB2_erf_int=quad(B2_erf_integrand, 0, np.inf, args=(KbT, pXi))
        pB2_diff=pB2_erf_int[0]-B2_WCA_int[0]
        #print(pB2_diff)
        if math.copysign(1, pB2_diff) > 0 :
            if sign==0 :
                B2_diff=pB2_diff
                B2_WCA_at_T=B2_WCA_int[0]
                B2_erf_at_T=pB2_erf_int[0]
                Xi_at_T=pXi 
                break 
        if math.copysign(1, pB2_diff) < 0 :
            if sign==1 :
                B2_diff=pB2_diff
                B2_WCA_at_T=B2_WCA_int[0]
                B2_erf_at_T=pB2_erf_int[0]
                Xi_at_T=pXi 
                break 
                #-----------------------------
                
                    
    #print "B2_diff=", B2_diff, "Xi_at_T=", Xi_at_T, "KbT=", KbT
    Xi.append(Xi_at_T)
    B2_erf.append(B2_erf_at_T)
    B2_WCA.append(B2_WCA_at_T)
    T.append(KbT)
    
    
##Plot B2_WCA vs T
plt.plot(T, B2_WCA, label = 'B2_WCA')
plt.plot(T, B2_erf, label = 'B2_erf')
plt.xlabel('KbT')
plt.ylabel('B2')
plt.title('B2 vs Temp')
plt.legend()

plt.show()


##Plot xi
plt.plot(T, Xi, label = 'Xi')
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
