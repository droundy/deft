#!/usr/bin/python3

#RUN this program from the directory it is listed in
#with command ./plot_WCA_Potential_fixed.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

#Plot WCA Potential vs r 
R=1/1.781797436  #for a sigma=1   DOESN"T WORK!!  graph wrong shape!
epsilon=1
sigma=1

#r=np.linspace(.4, .1, 2000)
#sigma_over_r_to_pow6=(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)*(sigma/r)
#V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon

#Plot repulsive portion from 0 to r=1.1225 (the min of the LJ potential):
r=np.linspace(.92, 1.11, 2000)  #repulsive portion only
r_dless=sigma/r  #plot dimensionless quantity!
sigma_over_r_to_pow6=(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)
V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon #LJ potential 
plt.plot(1/r_dless,V, label='LJ repulsive portion', color='blue')


#Plot repulsive portion from r=1.1225 to rmax:
r=np.linspace(1.11, 2.5, 2000) 
V=(r/r) - 1
plt.plot(r/sigma, V/epsilon, color='blue')


plt.xlabel('r/$\sigma$')
plt.ylabel('$V(r)/\epsilon$')
plt.title('WCA Potential')
#plt.legend()
plt.savefig('WCA_Potential.png')

# plt.show()



