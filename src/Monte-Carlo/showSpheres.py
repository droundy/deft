#!/usr/bin/python
from visual import *

allSpheres = genfromtxt("Spheres.dat")
spheres = allSpheres[0]
print(allSpheres)
print(spheres)

spheres = reshape(spheres, (-1, 3))
print(spheres)
for s in spheres:
    sphere(pos = s)
