#!/usr/bin/python2

import glob

f = open('figs/hughes-single-rod-in-water.dat', 'w')
for x in sorted(glob.glob('figs/hughes-single-rod-*nm-energy.dat')):
    with open(x) as r:
        f.write(r.read())

f = open('figs/single-rod-in-water.dat', 'w')
for x in sorted(glob.glob('figs/single-rod-*nm-energy.dat')):
    with open(x) as r:
        f.write(r.read())
