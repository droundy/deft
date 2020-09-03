#!/usr/bin/python3

#RUN this program from the directory it is listed in
#with command ./plot_LJ_Potential_2part.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

#Plot WCA Potential vs r 
#R=1/1.781797436  #for a sigma=1   DOESN'T WORK!!  graph wrong shape!
R=1/1.781797436 
epsilon=1
sigma=1


#Plot repulsive portion from 0 to r=1.1225 (the min of the LJ potential):
r=np.linspace(.92, 1.11, 2000)  #repulsive portion only
r_dless=sigma/r  #plot dimensionless quantity!
sigma_over_r_to_pow6=(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)
V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6) + epsilon #LJ potential 
plt.plot(1/r_dless,V, label='LJ repulsive portion', color='blue')

#Plot attractive portion from 0 to r=1.1225:
r=np.linspace(.92, 1.11, 2000) 
V=(r/r) - 2
plt.plot(1/r_dless,V, label='LJ attractive portion', color='red')

#Plot attractive portion from r=1.1225 to rmax:
r=np.linspace(1.11, 2.5, 2000) #attractive portion only
r_dless=sigma/r  #plot dimensionless quantity!
sigma_over_r_to_pow6=(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)*(r_dless)
V=4*epsilon*(sigma_over_r_to_pow6*sigma_over_r_to_pow6 - sigma_over_r_to_pow6)  #LJ potential 
plt.plot(1/r_dless,V, color='red')
plt.annotate('Attractive portion', xy=(1.45,-0.5), xytext=(1.6,-0.7),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)

#Plot repulsive portion from r=1.1225 to rmax:
r=np.linspace(1.11, 2.5, 2000) 
V=(r/r) - 1
plt.plot(1/r_dless,V, color='blue')
plt.annotate('Repulsive portion', xy=(1,2), xytext=(1.2,2.5),
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),)


plt.xlabel('r/$\sigma$')
plt.ylabel('$V(r)/\epsilon$')
plt.title('Leonard-Jones Potential in 2 parts')
#plt.legend()
plt.savefig("LJ_Potential_2part.pdf")

# plt.show()



