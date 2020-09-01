#!/usr/bin/python3

#RUN this program from the directory it is listed in
#with command ./plot_LJ_Potential.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

#Plot WCA Potential vs r 
#R=1/1.781797436  #for a sigma=1   DOESN'T WORK!!  graph wrong shape!
R=1/1.781797436 
epsilon=1
sigma=1

#print sigma

#r=np.linspace(.1, 2*R, 200)

#r=np.linspace(.9, 4, 200)  #SAVE!!! for plotting r
r=np.linspace(.9, 2.5, 200)  

r_dless=sigma/r  #plot dimensionless quantity!


sigma_over_r_to_pow6=(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)

#V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon #WCA potential
#V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6)  #LJ potential but looks like WCA
V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6)  #LJ potential 

plt.plot(1/r_dless,V)
plt.xlabel('r/$\sigma$')
#plt.xlabel('r')
plt.ylabel('V(r)/$\epsilon$')
plt.title('Leonard-Jones Potential')
#plt.legend()
plt.savefig("LJ_Potential.png")

# plt.show()



