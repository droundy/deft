#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy
import pylab

nm = 18.8972613

def plot_density(fin, col):
    gpermL=4.9388942e-3/0.996782051315 # conversion from atomic units to mass density
    x = []
    y = []
    for line in fin:
        current = str(line)
        pieces = current.split('\t')
        x.append(float(pieces[0])/nm)
        y.append(float(pieces[1])/gpermL)
    pyplot.plot(x, y, color = col)

#colors for plot
col1 = "#99ff88"
col2 = "#44dd55"
col3 = "#33bb88"
col4 = "#337799"
col5 = "#1133aa"
col6 = "#000044"

#make plot from files
ax1 = pyplot.subplot(2,1,1)
rod1 = open('figs/single-rod-slice-0.10.dat', 'r')
plot_density(rod1, col1)
rod2 = open('figs/single-rod-slice-0.30.dat', 'r')
plot_density(rod2, col2)
rod3 = open('figs/single-rod-slice-0.60.dat', 'r')
plot_density(rod3, col3)
rod4 = open('figs/single-rod-slice-1.00.dat', 'r')
plot_density(rod4, col4)
rod5 = open('figs/single-rod-slice-1.60.dat', 'r')
plot_density(rod5, col5)
rod6 = open('figs/single-rod-slice-2.00.dat', 'r')
plot_density(rod6, col6)
pyplot.hlines(1, 0, 1.3, 'k', '--')

ax2 = pyplot.subplot(2,1,2, sharex=ax1)
rod11 = open('figs/hughes-single-rod-slice-0.10.dat', 'r')
plot_density(rod11, col1)

rod12 = open('figs/hughes-single-rod-slice-0.30.dat', 'r')
plot_density(rod12, col2)
rod13 = open('figs/hughes-single-rod-slice-0.60.dat', 'r')
plot_density(rod13, col3)
rod14 = open('figs/hughes-single-rod-slice-1.00.dat', 'r')
plot_density(rod14, col4)
rod15 = open('figs/hughes-single-rod-slice-1.60.dat', 'r')
plot_density(rod15, col5)
rod16 = open('figs/hughes-single-rod-slice-2.00.dat', 'r')
plot_density(rod16, col6)
pyplot.hlines(1, 0, 1.3, 'k', '--')

data = pylab.loadtxt('figs/hughes-single-rod-slice-0.10.dat');
pylab.plot(data[:,0]/nm, data[:,3], 'r-')

data = pylab.loadtxt('figs/hughes-single-rod-slice-0.30.dat');
pylab.plot(data[:,0]/nm, data[:,3], 'm-')

data = pylab.loadtxt('figs/hughes-single-rod-slice-0.60.dat');
pylab.plot(data[:,0]/nm, data[:,3], 'b-')

#plot properties
pyplot.ylabel('                                      Density (g/mL)')
pyplot.xlabel('Radius (nm)')
pyplot.xlim(0, 1.3)
xticklabels = ax1.get_xticklabels()
pyplot.setp(xticklabels, visible=False)
pyplot.subplots_adjust(hspace=0.1)
pyplot.savefig('figs/density-compare.eps')
