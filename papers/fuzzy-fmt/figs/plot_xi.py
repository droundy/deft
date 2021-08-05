#!/usr/bin/python3

#This program plots Xi(T) and alpha(T)
#Run this program from deft/papers/fuzzy-fmt 
#with command ./figs/plot_xi.py

import numpy as np
import matplotlib.pyplot as plt

import find_xi

T=[]
Xi=[]
alpha=[]

for kT in np.arange(0.5, 200, 0.5):
    alpha_at_T = find_xi.find_alpha(kT)
    Xi_at_T = find_xi.find_Xi(alpha, kT)
    #print ('T=', kT)
    #print ('Xi_at_T=', Xi_at_T)
    #print ('alpha_at_T=', alpha_at_T)
    alpha.append(alpha_at_T)
    Xi.append(Xi_at_T)
    T.append(kT)
    
plt.figure("Xi(T)") 
    
##Plot Xi
plt.plot(T, Xi, label = 'Xi')
plt.xlabel('KbT')
plt.ylabel('Xi')
plt.title('Xi vs Temp')
plt.legend()

plt.figure("alpha(T)") 

##Plot alpha
plt.plot(T, alpha, label = 'alpha')
plt.xlabel('KbT')
plt.ylabel('alpha')
plt.title('alpha vs Temp')
plt.legend()

plt.show() 