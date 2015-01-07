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





if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)

allSpheres = genfromtxt(sys.argv[1])


balls = []
#for spheres in allSpheres:
spheres = allSpheres # reshape(spheres, (-1, 3))
#   rate(100)
bins=15

print 'next is lines'
divisions=numpy.zeros((bins,2))

for j in range(0,bins):
    divisions[j][0]=float(j/1.0)


'''
for k in range(1400,1450):
	print divisions[k][0]

x = spheres[1300][0]
print x
y = spheres[1300][1]
print y
z = spheres[1300][2]
print z
dist=pow((x*x+y*y+z*z),.5)
print round(dist,2)
if (dist<15):
	div=round(dist*100,0)
	divisions[div][1] =	divisions[div][1]+1
print divisions[div][1]


'''
print divisions[1:]
for i in range(len(spheres)):
	x = spheres[i][0]
	y = spheres[i][1]
	z = spheres[i][2]
	dist=pow((x*x+y*y+z*z),.5)
	round(dist,2)
	if (dist<14):
		div=round(dist*1.0)
		divisions[div][1] =	divisions[div][1]+1

#for k in range(1000,1450):
#	print divisions[k][0],divisions[k][1]
n=numpy.zeros((bins,1))
m=numpy.zeros((bins,1))
for l in range(0,bins):
	n[l]=divisions[l][0]
	m[l]=divisions[l][1]
	
lims = xlim()
xlim([1, 15]) 
plt.plot(n,m,'ro')
print n[10]

#savefig('test.pdf', bbox_inches=0)
plt.show()



