#!/usr/bin/env python

import math


fin = open('figs/single-rod-in-water.dat', 'r')
fout = open('figs/single-rods-calculated-density.dat', 'w')

kB = 3.16681539628059e-6 # This is Boltzmann's constant in Hartree/Kelvin
first = 1
nm = 18.8972613

for line in fin:
    current = str(line)
    pieces = current.split('\t')
    if first:
        r2 = float(pieces[0])/2*nm
        E2 = float(pieces[1])
        first = 0
    else:
        if ((float(pieces[0])/2*nm - r2) > 0.25):
            r1 = r2
            r2 = float(pieces[0])/2*nm
            E1 = E2
            E2 = float(pieces[1]) # actually it's energy per unit length!
            length = 1 # arbitrary
            r = (r1 + r2)/2
            dEdR = (E2-E1)/(r2-r1)*length
            area = 2*math.pi*r*length
            force = dEdR
            pressure = force/area
            kT = kB*298 # about this
            ncontact = pressure/kT
            fout.write(str(r)+'\t'+str(ncontact)+'\n')

fin.close()
fout.close()
