#!/usr/bin/env python

#need this to run without xserver
import matplotlib
matplotlib.use('Agg')

import math
import matplotlib.pyplot as pyplot
import numpy

def plot_X(fin, col, lin):
    nm = 18.8972613
    x = []
    y = []
    for line in fin:
        current = str(line)
        pieces = current.split('\t')
        x.append(float(pieces[0])/nm)
        y.append(4.0*(1-float(pieces[2])))
    pyplot.plot(x, y, color = col, linestyle = lin)

#colors for plot
col1 = "#99ff88"
col2 = "#44dd55"
col3 = "#33bb88"
col4 = "#337799"
col5 = "#1133aa"
col6 = "#000044"

oldline = "-"
newline = ":"

#make plot from files
rod1 = open('figs/single-rod-slice-00.1.dat', 'r')
plot_X(rod1, col1, newline)
rod2 = open('figs/single-rod-slice-00.3.dat', 'r')
plot_X(rod2, col2, newline)
rod3 = open('figs/single-rod-slice-00.6.dat', 'r')
plot_X(rod3, col3, newline)
rod4 = open('figs/single-rod-slice-01.0.dat', 'r')
plot_X(rod4, col4, newline)
rod5 = open('figs/single-rod-slice-01.6.dat', 'r')
plot_X(rod5, col5, newline)
rod6 = open('figs/single-rod-slice-02.0.dat', 'r')
plot_X(rod6, col6, newline)

rod11 = open('figs/hughes-single-rod-slice-00.1.dat', 'r')
plot_X(rod11, col1, oldline)
rod12 = open('figs/hughes-single-rod-slice-00.3.dat', 'r')
plot_X(rod12, col2, oldline)
rod13 = open('figs/hughes-single-rod-slice-00.6.dat', 'r')
plot_X(rod13, col3, oldline)
rod14 = open('figs/hughes-single-rod-slice-01.0.dat', 'r')
plot_X(rod14, col4, oldline)
rod15 = open('figs/hughes-single-rod-slice-01.6.dat', 'r')
plot_X(rod15, col5, oldline)
rod16 = open('figs/hughes-single-rod-slice-02.0.dat', 'r')
plot_X(rod16, col6, oldline)

#plot properties
pyplot.ylabel('Number of hydrogen bonds')
pyplot.xlabel('Radius (nm)')
pyplot.xlim(0, 1.3)
pyplot.savefig('figs/single-rod-X-plot.eps')
