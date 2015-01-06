#!/usr/bin/python
from visual import *
import sys
import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib.widgets import Slider, Button, RadioButtons
sys.setrecursionlimit(2500)


def fart(str):
	ax.cla()
	x0=0
	y0=0
	z0=0
	ax.axhline(y=5,xmin=-10,xmax=10,linewidth=1, color='r')
	ax.axhline(y=5.5,xmin=-10,xmax=10,linewidth=1, color='r')
	for i in range(len(spheres)):
		x = spheres[i][0]
		y = spheres[i][1]
		z = spheres[i][2]
	        if (z < scoord.val and z > (scoord.val -.5)):
		    r = (z+1)-scoord.val
		    if (y<5.5 and y>5):
			cir = plt.Circle((x,y), radius=r,  fc='b')

	            else: cir = plt.Circle((x,y), radius=r,  fc='y')
		    ax.add_patch(cir)
	if (z0 < scoord.val and z0 > (scoord.val -.5)):
	    print 'pofpdofpdakf'
	    r = (z0+1)-scoord.val
	    cir2 = plt.Circle((x0,y0), radius=r,  fc='r')
	    ax.add_patch(cir2)

	


if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)

allSpheres = genfromtxt(sys.argv[1])


balls = []
#for spheres in allSpheres:
spheres = allSpheres # reshape(spheres, (-1, 3))
#   rate(100)

ax = subplot(111, autoscale_on=False, aspect='equal') 
ax.set_xbound(-20,20)
ax.set_ybound(-20,20)
subplots_adjust(left=0.25, bottom=0.25)

a0 = 0
axcoord = axes([0.25,.01,.65,.03])
scoord = Slider(axcoord, 'x cord', -20, 20, valinit=a0)
print 'next is lines'
print spheres[3][1]
print str(scoord.val)

for i in range(len(spheres)):
	x = spheres[i][0]
	y = spheres[i][1]
	z = spheres[i][2]
        if (z < scoord.val and z > (scoord.val -.01)):
            r = (z+1)+10*scoord.val
	    cir = plt.Circle((x,y), radius=r,  fc='y')
	    ax.add_patch(cir)
	    
  #      elif (c == 0):
  #          cir = plt.Circle((x,y), radius=r,  fc='r')
        #else:
   #     cir = plt.Circle((x,y), radius=r, color='b', fill=False)
        
scoord.on_changed(fart)
  
plt.axhline(y=1,xmin=0,xmax=10,color='r')
plt.axhline(y=1.2,xmin=0,xmax=10,color='r')
savefig('test.pdf', bbox_inches=0)
plt.show()



