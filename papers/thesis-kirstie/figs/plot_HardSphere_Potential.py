#!/usr/bin/python3

#RUN this program from the directory it is listed in
#with command ./plot_HardSphere_Potential.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

#Plot Hard Sphere Potential vs r 
R=1/1.781797436  #for a sigma=1   DOESN"T WORK!!  graph wrong shape!
epsilon=1
sigma=1

#Plot zero portion from 0 to r=2R=1.1225:  #to shift plot to right
r=np.linspace(0.01, 1.22, 2000) #zero only
V=(r/r)-1 
plt.plot(r/sigma, V/epsilon, color='white')

#Plot infinte portion from 0 to r=2R=1.1225 :
V=np.linspace(0.01, 3, 2000)  #infinite portion only
r=(V/V)+0.1225 
plt.plot(r/sigma,V/epsilon, color='blue')

#Plot zero portion from r=2R=1.1225 to rmax:
r=np.linspace(1.1225, 4, 2000) #zero portion only
V=(r/r)-1   #the r/r makes x and y have the same dimensions for plotting
plt.plot(r/sigma, V/epsilon, color='blue')


#plt.xlabel('r/$\sigma$')
#plt.ylabel('V_HS/$\epsilon$')
plt.xlabel('r')
#plt.ylabel('$V_{HS}(r)/\epsilon$')
plt.ylabel('V(r)')
plt.title('Hard Sphere Potential')
#plt.legend()
plt.savefig("plot_HardSphere_Potential.pdf")

# plt.show()



