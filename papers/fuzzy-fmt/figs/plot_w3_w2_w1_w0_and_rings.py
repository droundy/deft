#!/usr/bin/python3

#RUN this program from the directory it is listed in
#with command ./plot_w3_w2_w1_w0_and_rings.py

#Note: The weight functions w2, w1, and w0 are spherical shells (or rings in 2D)
#with max value at alpha/2 for w2 (when centered about origin) and slightly below
#alpha/2 for w1 and even more slightly below alpha/2 for w0. 
#Weight function w3 is a sphere of constant value that tapers off at the edge, 
#and the middle of the tapering occurs at alpha/2.

##Weight functions are calculated for rprime=0, and magnitude of r=rz (where rx=ry=0).

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

import find_xi

sigma=1
epsilon=1
PI=np.pi

def w3(r):
    return (1-special.erf((r-(alpha/2))/Xi))/2  

def w2(r):
    return (1/(Xi*np.sqrt(PI)))*np.exp(-(((r-(alpha/2))/Xi)*((r-(alpha/2))/Xi)))  
    
def w1(r):
    return w2(r)/(4*PI*r)

def w0(r):
    return w2(r)/(4*PI*r*r)


plt.figure('weights vs r')

##Plot w0, w1, w2, w3 vs rz at KbT
KbT=2    #good up to 200
alpha= find_xi.find_alpha(KbT)
Xi = find_xi.find_Xi(alpha, KbT)
rz=np.linspace(.01, 1, 20000)
plt.plot(rz/(alpha/2),w0(rz), label = 'w0', color='purple')
plt.plot(rz/(alpha/2),w1(rz), label = 'w1', color='red')
plt.plot(rz/(alpha/2),w2(rz), label = 'w2', color='blue')
plt.plot(rz/(alpha/2),w3(rz), label = 'w3', color='green')
plt.xlim(0, 2)
plt.xlabel('rz/(alpha/2)')
plt.ylabel('weight function')
plt.title('Weight Functions at KbT=%g' % KbT)
plt.legend()
plt.savefig('weight_functions.pdf')


plt.figure('w2 ring')

##Plot w2 ring at KbT
x = np.arange(-1,1,0.01)
X,Y = np.meshgrid(x,x)
plt.pcolormesh(X,Y, w2(np.sqrt(X**2+Y**2)))
plt.colorbar()
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('w2-ring.pdf')


# plt.figure('w3 ring')  #doesn't work

##Plot w3 shpere at KbT
# x = np.arange(-1,1,0.01)
# X,Y = np.meshgrid(x,x)
# plt.pcolormesh(X,Y, w3(np.sqrt(X**2+Y**2)))
# plt.colorbar()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.savefig('w3-ring.pdf')

plt.show()
