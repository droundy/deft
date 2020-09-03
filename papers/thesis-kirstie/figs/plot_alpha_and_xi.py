#!/usr/bin/python3

#This program plots alpha derived from matching Verf and Vwca at r=alpha
#and plots xi derived from matching dVerf/dr and dVwca/dr at r=alpha

#RUN this program from the directory it is listed in
#with command ./plot_alpha_and_xi.py


from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math


sigma=1
epsilon=1
PI=np.pi


#KbT=np.linspace(.05, 200, 2000) #transistion
KbT=np.linspace(.05, 2000, 2000) 
#KbT=np.linspace(.05, 1000000, 2000) 
#KbT=np.linspace(.05, 1e50, 2000) 
 
alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))  
xi=alpha/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/KbT)+np.log(2)))

##Plot alpha
plt.plot(KbT,alpha, label = 'alpha')
plt.xlabel('KbT')
plt.ylabel('alpha')
plt.title('Alpha vs Temp')
#plt.legend()
plt.savefig('plot_alpha.pdf')

plt.figure()

##Plot xi
plt.plot(KbT,xi, label = 'xi')
plt.xlabel('KbT')
plt.ylabel('xi')
plt.title('Xi vs Temp')
#plt.legend()
plt.savefig('plot_xi.pdf')

#plt.show()
