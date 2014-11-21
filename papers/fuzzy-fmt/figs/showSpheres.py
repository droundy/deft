#!/usr/bin/python
from visual import *
import sys

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat")
    exit(1)

allSpheres = genfromtxt(sys.argv[1])
#print allSpheres.size/3
#print(allSpheres)
#print(spheres)

#container = box(opacity=.1, color=color.red)
#container.size=(6,6,6)
#container = sphere(opacity=.1, color=color.red)
#container.radius=(3)

balls = []
#for spheres in allSpheres:
spheres = allSpheres # reshape(spheres, (-1, 3))
#   rate(100)
for i in range(len(spheres)):
   if i >= len(balls):
      balls.append(sphere(pos=spheres[i]))
   else:
      balls[i].pos = spheres[i]
   for j in range(len(spheres)):
     if i != j:
         d = sqrt((spheres[j][0]-spheres[i][0])**2+(spheres[j][1]-spheres[i][1])**2+(spheres[j][2]-spheres[i][2])**2)
         if d < 2:
           print 'distance is', i, j, d
           balls[i].color = color.yellow

#The rest is stuff that Jeff added for investigating a bug with the mc-walls data.  More info about walls!
countOutside = 0
countInside = 0
countOutsideZ = 0
spheresOutside = []
for i in range (len(spheres)):
    if ((abs(spheres[i][0])<10) & (abs(spheres[i][1])<10) & (abs(spheres[i][2])<10)):
        #print(spheres[i])
        countInside += 1
    else:
        if (abs(spheres[i][2])>10):
            countOutsideZ += 1
        spheresOutside.append(spheres[i])
        countOutside += 1

print ("And the spheres outside are: ")
#print (spheresOutside)
print ("count inside: " + str(countInside) + "  And count outside: " + str(countOutside))
print ("And the number that are outside in Z coordinate are: " + str(countOutsideZ))        

