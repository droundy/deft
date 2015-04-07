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





#if len(sys.argv) < 2:
#    print("Usage:  " + sys.argv[0] + " filename.dat")
#    exit(1)

#allSpheres = genfromtxt(sys.argv[1])


balls = []
#for spheres in allSpheres:
#spheres = allSpheres # reshape(spheres, (-1, 3))
#   rate(100)

ax = subplot(111, autoscale_on=False, aspect='equal') 
ax.set_xbound(-20,20)
ax.set_ybound(-20,20)
subplots_adjust(left=0.25, bottom=0.25)

a0 = 0
#axcoord = axes([0.25,.01,.65,.03])


r = 1
cir = plt.Circle((0,0), radius=r,  fc='b')
ax.add_patch(cir)
cir = plt.Circle((0,0), radius=5,  color='r',fill=False)
ax.add_patch(cir)
cir = plt.Circle((0,0), radius=6,  color='r',fill=False)
ax.add_patch(cir)
colr='y'
for i in range(0,35):
	x = randint(-20,20)
	y = randint(-20,20)
	if pow(x*x + y*y,.5) < 6 and pow(x*x + y*y,.5) > 5:
		colr = 'r'
	cir = plt.Circle((x,y), radius=r,  fc=colr)
	ax.add_patch(cir)
	colr='y'
   

  

savefig('test2.pdf', bbox_inches=0)
plt.show()



