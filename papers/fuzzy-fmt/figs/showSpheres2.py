#!/usr/bin/python
import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
sys.setrecursionlimit(2500)
if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)


ax = subplot(111, autoscale_on=False, aspect='equal') 
ax.set_xbound(-25,25)
ax.set_ybound(-25,25)
subplots_adjust(left=0.25, bottom=0.25)

a0 = 0
axcoord = axes([0.25,.1,.65,.03])
scoord = Slider(axcoord, 'x cord', -20, 20, valinit=a0)


allSpheres = numpy.loadtxt(sys.argv[1])




balls = []

spheres = allSpheres # reshape(spheres, (-1, 3))



    
lx = [[None]*4 for m in range(len(spheres))]
inter = 0
lenz = 20
j = a0

#for l in range (len(spheres)):
#    for j in range(len(spheres)):
#        if l != j:
#            d = sqrt((spheres[j][0]-spheres[l][0])**2+(spheres[j][1]-spheres[l][1])**2+(spheres[j][2]-spheres[l][2])**2)
#            if d < 2:
#                print 'distance is ', d
#                lx[j][3]=1
#print 'done'
            
                

def plotDat(j):
    ax.clear()
    inter = 0
    for k in range (len(spheres)):   
        z = spheres[k][2]
        if (abs(j-z) < 1):
            # print spheres[k][0], spheres[k][1], spheres[k][2]
            x = (spheres[k][0])
            y = (spheres[k][1])
            z = (spheres[k][2])
            lx[inter][0] = x
            lx[inter][1] = y
            si = abs(j-z)
            lx[inter][2] = sqrt(1-(si*si)) 
            inter = inter +1

#        print 'in z ', k, ' ', len(spheres)
    print 'inter is ' ,inter
    for i in range(inter):
        x = lx[i][0]
        y = lx[i][1]
        r = lx[i][2]
       # c = lx[i][3]
        if (lx[i][3] == 1):
            cir = plt.Circle((x,y), radius=r,  fc='y')
  #      elif (c == 0):
  #          cir = plt.Circle((x,y), radius=r,  fc='r')
        else:
            cir = plt.Circle((x,y), radius=r, fc = 'b')
        ax.add_patch(cir)    

def update(val):    
    ax.clear()
    plotDat(scoord.val)
    show()
    
scoord.on_changed(update)
plotDat(scoord.val)
update(a0)






