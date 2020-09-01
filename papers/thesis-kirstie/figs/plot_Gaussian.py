#!/usr/bin/python2

#RUN this program from the directory it is listed in
#with command ./plot_Gaussian.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

plot1="plot_Gauss.png"

#Plot Gaussian vs rprime_magnitude at fixed gw  R=0
rprime=np.linspace(0, .01, 20000)
gw=.001
 
#Gaussian=(1/np.sqrt(2*np.pi)*gw)*(1/np.sqrt(2*np.pi)*gw)*(1/np.sqrt(2*np.pi)*gw)*np.exp(-(rprime/np.sqrt(2)*gw)*(rprime/np.sqrt(2)*gw))

Gaussian=np.exp(-(rprime/(np.sqrt(2)*gw))*(rprime/(np.sqrt(2)*gw)))

plt.plot(rprime,Gaussian)
plt.xlabel('rprime')
plt.ylabel('Gaussian')
plt.title('Gaussian at gw=0.001')
#plt.legend()
plt.savefig(plot1)

plt.show()



