#!/usr/bin/env python

from __future__ import division
import matplotlib.pyplot as pyplot
import numpy, math

#Get data before transition
lines = 0
f = open('rods-slice-01.0-01.1.dat','r')
for line in f:
    lines +=1
f.close()

rowcol = math.sqrt(lines)
x = numpy.zeros((rowcol, rowcol))
y = numpy.zeros((rowcol, rowcol))
density_1 = numpy.zeros((rowcol, rowcol))
row = 0
col = 0

f = open('rods-slice-01.0-01.1.dat','r')
for line in f:
    current = str(line)
    pieces = current.split ('\t')
    if col == (rowcol):
        col = 0
        row +=1

    x[row][col] = float(pieces[1])
    y[row][col] = float(pieces[2])
    density_1[row][col] = float(pieces[3])
    col +=1

f.close()

#get data for after transition
lines = 0
f = open('rods-slice-01.0-01.4.dat','r')
for line in f:
    lines +=1
f.close()

rowcol = math.sqrt(lines)
a = numpy.zeros((rowcol, rowcol))
b = numpy.zeros((rowcol, rowcol))
density_2 = numpy.zeros((rowcol, rowcol))
row = 0
col = 0

f = open('rods-slice-01.0-01.4.dat','r')
for line in f:
    current = str(line)
    pieces = current.split ('\t')
    if col == (rowcol):
        col = 0
        row +=1

    a[row][col] = float(pieces[1])
    b[row][col] = float(pieces[2])
    density_2[row][col] = float(pieces[3])
    col +=1
f.close()

pyplot.set_cmap('GnBu')

pyplot.subplot(211)
pyplot.imshow(density_1)
#pyplot.pcolor(x,y,density_1)

pyplot.subplot(212)
pyplot.imshow(density_2)
#pyplot.pcolor(a,b,density_2)

pyplot.subplots_adjust(bottom=0.1, right=1.0, top=0.9, left=0.3)
cax = pyplot.axes([0.85, 0.1, 0.045, 0.8])
pyplot.colorbar(cax=cax)
pyplot.show()
