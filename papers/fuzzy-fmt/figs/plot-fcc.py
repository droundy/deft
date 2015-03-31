#!/usr/bin/python
from visual import *
import sys
import pylab, numpy, sys, math
from pylab import *
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 2:
    print("Usage:  " + sys.argv[0] + " filename.dat.fcc")
    exit(1)


ax = subplot(111, autoscale_on=False, aspect='equal') 
ax.set_xbound(-25,25)
ax.set_ybound(-25,25)
subplots_adjust(left=0.25, bottom=0.25)


data = genfromtxt(sys.argv[1])


def plotData():
    ax.clear
    for i in range(len(data)):
        if data[2] != 0:
            cir=plt.Circle((data[0],data[1]),radius=1,fc='y')
            ax.add_patch(cir)
show()


