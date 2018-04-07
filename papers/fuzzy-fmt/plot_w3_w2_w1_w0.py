#!/usr/bin/python2

#RUN this program from the directory it is listed in
#with command ./plot_w3_w2_w1_w0.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math


plot1="plot_w0_w1_w2_w3_vs_KbT.png"
plot2="plot_w0_w1_w2_w3_vs_rz.png"

sigma=1
epsilon=1
PI=np.pi

#Note: the weight functions w2, w1, and w0 are spherical shells (or rings in 2D)
#with max value at alpha/2 for w2 (when centered about origin) and slightly below
#alpha/2 for w1 and even more slightly below alpha/2 for w0. 
#Weight function w3 is a sphere of constant value that tapers off at the edge, 
#and the middle of the tapering occurs at alpha/2.

#Weight functions are calculated for rprime=0, and magnitude of r=rz (where rx=ry=0).

#Plot w0, w1, w2, w3 vs KbT at rz=.55 
KbT=np.linspace(.0005, 1, 20000) #replace with a range of values
rz=.55
 
alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))  
zeta=alpha/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/KbT)+np.log(2)))
w2=(1/(zeta*np.sqrt(PI)))*np.exp(-(((rz-(alpha/2))/zeta)*((rz-(alpha/2))/zeta)))
w1=w2/(4*PI*rz)
w0=w2/(4*PI*rz*rz)
w3=(1-special.erf((rz-(alpha/2))/zeta))/2 

plt.plot(KbT,w2, label = 'w2')
plt.plot(KbT,w1, label = 'w1')
plt.plot(KbT,w0, label = 'w0')
plt.plot(KbT,w3, label = 'w3')
plt.xlabel('KbT')
plt.ylabel('weight function')
plt.title('Weight Functions vs Temp at rz=.55')
plt.legend()
plt.savefig(plot1)

plt.show()


#Plot w0, w1, w2, w3 vs rz at KbT=2
KbT=2
rz=np.linspace(.01, 1, 20000)

alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))  
zeta=alpha/(6*np.sqrt(PI)*(np.sqrt((epsilon*np.log(2))/KbT)+np.log(2)))
w2=(1/(zeta*np.sqrt(PI)))*np.exp(-(((rz-(alpha/2))/zeta)*((rz-(alpha/2))/zeta)))
w1=w2/(4*PI*rz)
w0=w2/(4*PI*rz*rz)
w3=(1-special.erf((rz-(alpha/2))/zeta))/2 

plt.plot(rz,w0, label = 'w0')
plt.plot(rz,w1, label = 'w1')
plt.plot(rz,w2, label = 'w2')
plt.plot(rz,w3, label = 'w3')
plt.xlabel('rz')
plt.ylabel('weight function')
plt.title('Weight Functions vs rz at KbT=2')
plt.legend()
plt.savefig(plot2)

plt.show()
