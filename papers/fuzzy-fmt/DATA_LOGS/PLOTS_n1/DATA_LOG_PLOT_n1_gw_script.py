#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt with command
# ./DATA_LOG_PLOT_script.py [directory where data stored] --f[fixed quantity] [value of fixed quantity] --x[] --y[]
#to create plots from plot.dat files already in the data directory
####ie. ENTER ./DATA_LOGS/PLOTS_n1/DATA_LOG_PLOT_script.py  pears --ycphi --xdx  --ptname addsomethingtoplotname   
#ie. ENTER ./DATA_LOGS/PLOTS_n1/DATA_LOG_PLOT_script.py
#to plot  diff_free_enery vs gw at fixed T=2

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

#parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")
#groupf = parser.add_mutually_exclusive_group(required=True)
#groupx = parser.add_mutually_exclusive_group(required=True)
#groupy = parser.add_mutually_exclusive_group(required=True)

 ####parser.add_argument('directory', metavar='exisiting_data_directory', type=str,
 ####                    help='exisiting directory for data files') 
 ####parser.add_argument('data_file', metavar='file with plot data', type=str,
 ####                    help='file with plot data') 
 ####groupf.add_argument('--ftemp', action="store_true",
 ####                    help='use plot.dat file with this fixed temperature')
 ####groupf.add_argument('--fdensity', action="store_true",
 ####                    help='use plot.dat file with this fixed temperature') 
 ####parser.add_argument('value', metavar='value_of_fixed_quantity', type=str,
 ####                    help='use plot.dat file with this fixed value')


#parser.add_argument('--xlab', metavar='label for x-axis', type=str,
#                    help='label for x-axis. use with --xcol.') 
#parser.add_argument('--ylab', metavar='label for y-axis', type=str,
#                    help='label for y-axis. use with --ycol.') 
 ####parser.add_argument('--ptname', metavar='include in name of plot', type=str,
 ####                    help='info added to plot name') 
                    
                
#groupx.add_argument('--xdx', action="store_true",
#                    help='dx on x-axis') 
#groupx.add_argument('--xcphi_1', action="store_true",
#                    help='cphi_1 on x-axis') 
#groupx.add_argument('--xcphi_2', action="store_true",
#                    help='cphi_2 on x-axis') 
#groupx.add_argument('--xcphi_3', action="store_true",
#                    help='cphi_3 on x-axis') 
#groupx.add_argument('--xcFE_vol', action="store_true",
#                    help='Total crystal Free Energy per volume on x-axis')   
#groupx.add_argument('--xcol', metavar='column for x-axis data', type=int,
#                    help='column for x-axis data') 

#groupy.add_argument('--ydx', action="store_true",
#                    help='dx on y-axis') 
#groupy.add_argument('--ycphi_1', action="store_true",
#                    help='cphi_1 on y-axis') 
#groupy.add_argument('--ycphi_2', action="store_true",
#                    help='cphi_2 on y-axis') 
#groupy.add_argument('--ycphi_3', action="store_true",
#                    help='cphi_3 on y-axis') 
#groupy.add_argument('--ycFE_vol', action="store_true",
#                    help='Total crystal Free Energy per volume on y-axis')   
#groupy.add_argument('--ycol', metavar='column for y-axis data', type=int,
#                    help='column for y-axis data') 

#args=parser.parse_args()
#print
#print args

####data_directory=args.directory

####if args.ftemp:
####    fixed_quantity="kT"
####elif args.fdensity:
####    fixed_quantity="n"
    
####fixed_value=args.value

#data_file=data_directory+"/plot_"+fixed_quantity+fixed_value+".dat"
#data_file=data_directory+"/plot.dat"
#thisdata = np.loadtxt(data_file)

#n=args.n

thisdata_n1_80000_gw = np.loadtxt("DATA_LOGS/PLOTS_n1/DATA_LOG_PLOT_n1_gw.dat")

MC80000_gw=1

if MC80000_gw > 0 :
  plt.figure()
  plt.plot(thisdata_n1_80000_gw[:,0], thisdata_n1_80000_gw[:,1], label = 'FEdiff')
  plot_name="DATA_LOGS/PLOTS_n1/plot_DATA_LOG_MC80000_gw"
  plot_title="FEdiff vs gw"
  plt.title(plot_title)
  plt.xlabel("gw")
  plt.ylabel("FEdiff")
  plt.legend()
  plt.savefig(plot_name)

plt.show()

