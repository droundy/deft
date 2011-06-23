#!/usr/bin/python
from visual import *

allSpheres = genfromtxt("Spheres.dat")
print allSpheres.size
print(allSpheres)
#print(spheres)

container = box(opacity=.1, color=color.red)
container.size=(5,5,5)

for spheres in allSpheres:
   spheres = reshape(spheres, (-1, 3))
   for s in spheres:
        sphere(pos = s)
   print(spheres)
