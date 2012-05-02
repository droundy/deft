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

