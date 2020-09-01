#!/usr/bin/python2

#RUN this program from the directory it is listed in
#with command ./plot_Gaussian.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

#Plot Gaussian vs rprime_magnitude at fixed gw  R=0
gw=1
x=np.linspace(-4*gw, 4*gw, 20000)
 
#Gaussian=(1/np.sqrt(2*np.pi)*gw)*(1/np.sqrt(2*np.pi)*gw)*(1/np.sqrt(2*np.pi)*gw)*np.exp(-(rprime/np.sqrt(2)*gw)*(rprime/np.sqrt(2)*gw))

def Gaussian(x):
    return np.exp(-(x/(np.sqrt(2)*gw))*(x/(np.sqrt(2)*gw)))

scale = 0.9
plt.figure(figsize=(4*scale,3*scale))

plt.xlabel('$x$')
plt.xticks([-gw,0,gw], ['$-\sigma$', '$0$', '$\sigma$'])

plt.plot([-gw, -gw], [0, Gaussian(gw)], 'k-')
plt.plot([gw, gw], [0, Gaussian(gw)], 'k-')

plt.ylim(0, 1.1)
plt.xlim(-3*gw, 3*gw)

plt.axvline(0, color='k')

plt.plot(x,Gaussian(x))

plt.yticks([0, 1], ['$0$', r'$\frac{1}{\sqrt{2\pi}\sigma}$'])
plt.axhline(1, linestyle=':', color='k')

# plt.annotate('hello', xy=(0,1), xytext=(0.5*gw,1),
#              arrowprops=dict(arrowstyle="->",
#                             connectionstyle="arc3"),)
plt.text(-2.5*gw, 0.75, r"$\frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{x^2}{2\sigma^2}}$")

plt.tight_layout()
#plt.legend()
plt.savefig("Gaussian.pdf")

# plt.show()



