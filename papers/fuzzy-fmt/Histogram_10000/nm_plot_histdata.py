#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./nm_plot_histdata.py [filename.dat]


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")
parser.add_argument('filedat', metavar='datafile', type=str,
                    help='file with data to plot') 
args=parser.parse_args()

thisdata = np.loadtxt(args.filedat)
#thisdata = np.loadtxt('Hist_10000.dat')

num_lines=25

for i in range(0, num_lines) :
   if (thisdata[i,1] == .01) :
        x =thisdata[i,2]   #selects mc_constant
        y =thisdata[i,0]   #selects percent_negative
        print y, x, '.01'
        plt.scatter(x,y, color='blue', label='gw=.01')
        
for i in range(0, num_lines) :
   if (thisdata[i,1] == .05) :
        x =thisdata[i,2]   #selects mc_constant
        y =thisdata[i,0]   #selects percent_negative
        print y, x, '.1'
        plt.scatter(x,y, color='lawngreen', label='gw=.05')

for i in range(0, num_lines) :
   if (thisdata[i,1] == .1) :
        x =thisdata[i,2]   #selects mc_constant
        y =thisdata[i,0]   #selects percent_negative
        print y, x, '.1'
        plt.scatter(x,y, color='red', label='gw=.1')
        
for i in range(0, num_lines) :
   if (thisdata[i,1] == .5) :
        x =thisdata[i,2]   #selects mc_constant
        y =thisdata[i,0]   #selects percent_negative
        print y, x, '.1'
        plt.scatter(x,y, color='magenta', label='gw=.5')

      
plt.axhspan(60, 40, color='green', alpha=0.15, lw=0)
plt.axhspan(49.5, 50.5, color='black', alpha=0.15, lw=0)

plt.title('Histogram results: gw0.01=blue, gw0.05=green, gw0.1=red, gw.5=purple')
plt.ylabel('percent of negative relerrors') 
plt.xlabel('mc_constant')       
#plt.xlabel('mc_constant (with mc_prefactor=10,000)')
plt.show()
