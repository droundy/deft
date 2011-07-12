#!/usr/bin/python
from visual import *

allSpheres = genfromtxt("Spheres.dat")
print allSpheres.size/3
print(allSpheres)
#print(spheres)

#container = box(opacity=.1, color=color.red)
#container.size=(6,6,6)
container = sphere(opacity=.1, color=color.red)
container.radius=(3)
for spheres in allSpheres:
   balls = []
   spheres = reshape(spheres, (-1, 3))
   for s in spheres:
         balls.append(sphere(pos = s))
   rate(1)
   for i in balls:
      i.visible = False
   del balls
   print(spheres)
