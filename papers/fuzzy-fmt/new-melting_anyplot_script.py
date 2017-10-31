#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt with command
# ./new-melting_anyplot_script.py [directory where data stored] --f[] [temp] --x[] --y[]
#to create plots from plot.dat files already in the data directory

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot.', epilog="stuff...")
groupf = parser.add_mutually_exclusive_group(required=True)
groupx = parser.add_mutually_exclusive_group(required=True)
groupy = parser.add_mutually_exclusive_group(required=True)

parser.add_argument('directory', metavar='exisiting_data_directory', type=str,
                    help='exisiting directory for data files') 
groupf.add_argument('--ftemp', action="store_true",
                    help='use plot.dat file with this fixed temperature')
groupf.add_argument('--fdensity', action="store_true",
                    help='use plot.dat file with this fixed temperature') 
parser.add_argument('value', metavar='value_of_fixed_quantity', type=str,
                    help='use plot.dat file with this fixed value')
parser.add_argument('--xlab', metavar='label for x-axis', type=str,
                    help='label for x-axis. use with --xcol.') 
parser.add_argument('--ylab', metavar='label for y-axis', type=str,
                    help='label for y-axis. use with --ycol.') 
                    
                     
groupx.add_argument('--xtemp', action="store_true",
                    help='temperature on x-axis') 
groupx.add_argument('--xdensity', action="store_true",
                    help='density on x-axis') 
groupx.add_argument('--xcfe', action="store_true",
                    help='crystal free energy/atom on x-axis') 
groupx.add_argument('--xhfe', action="store_true",
                    help='homogeneous free energy/atom on x-axis') 
groupx.add_argument('--xdiff', action="store_true",
                    help='diff in free energy on x-axis')   
groupx.add_argument('--xcfev', action="store_true",
                    help='crystal energy/volume on x-axis') 
groupx.add_argument('--xfv', action="store_true",
                    help='fraction of vacancies (fv) on x-axis') 
groupx.add_argument('--xgw', action="store_true",
                    help='Gaussian width on x-axis')
groupx.add_argument('--xcol', metavar='column for x-axis data', type=int,
                    help='column for x-axis data') 

groupy.add_argument('--ytemp', action="store_true",
                    help='temperature on y-axis') 
groupy.add_argument('--ydensity', action="store_true",
                    help='density on y-axis') 
groupy.add_argument('--ycfe', action="store_true",
                    help='crystal free energy/atom on y-axis') 
groupy.add_argument('--yhfe', action="store_true",
                    help='homogeneous free energy/atom on y-axis') 
groupy.add_argument('--ydiff', action="store_true",
                    help='diff in free energy on y-axis')   
groupy.add_argument('--ycfev', action="store_true",
                    help='crystal energy/volume on y-axis') 
groupy.add_argument('--yfv', action="store_true",
                    help='fraction of vacancies (fv) on y-axis') 
groupy.add_argument('--ygw', action="store_true",
                    help='Gaussian width on y-axis')
groupy.add_argument('--ycol', metavar='column for y-axis data', type=int,
                    help='column for y-axis data') 

args=parser.parse_args()
#print
#print args

data_directory=args.directory

if args.ftemp:
    fixed_quantity="kT"
elif args.fdensity:
    fixed_quantity="rd"
    
fixed_value=args.value

data_file=data_directory+"/plot_"+fixed_quantity+fixed_value+".dat"
thisdata = np.loadtxt(data_file)
print
print "Using data from file:"+data_file
print

if args.xtemp:  
    x_axis=thisdata[:,0]    
    x_label="Temperature (kT)"
    x_plot="kT"
elif args.xdensity:  
    x_axis=thisdata[:,1]     
    x_label="Reduced Density (n)"
    x_plot="n"
elif args.xcfe:  
    x_axis=thisdata[:,2]     
    x_label="Crystal Free Energy/atom"
    x_plot="cFE"
elif args.xhfe:  
    x_axis=thisdata[:,3]     
    x_label="Homogeneous Free Energy/atom" 
    x_plot="hFE"
elif args.xdiff:  
    x_axis=thisdata[:,4]     
    x_label="Diff=(CrystalFE-HomogeneousFE)/atom" 
    x_plot="DiffFE"
elif args.xcfev:  
    x_axis=thisdata[:,5]     
    x_label="Crystal Free Energy/volume"
    x_plot="cFEv"
elif args.xfv:  
    x_axis=thisdata[:,6]     
    x_label="Fraction of vacancies (fv)"
    x_plot="fv"
elif args.xgw:   
    x_axis=thisdata[:,7]     
    x_label="Width of Gaussian (gwidth)"
    x_plot="gw"
elif args.xcol:   
    x_axis=thisdata[:,args.xcol]     
    x_label=args.xlab
    x_plot=args.xlab
    
if args.ytemp:  
    y_axis=thisdata[:,0]     
    y_label="Temperature (kT)"
    y_plot="kT"
elif args.ydensity:  
    y_axis=thisdata[:,1]     
    y_label="Reduced Density (n)"
    y_plot="n"
elif args.ycfe:  
    y_axis=thisdata[:,2]     
    y_label="Crystal Free Energy/atom"
    y_plot="cFE"
elif args.yhfe:  
    y_axis=thisdata[:,3]    
    y_label="Homogeneous Free Energy/atom"
    y_plot="hFE" 
elif args.ydiff:  
    y_axis=thisdata[:,4]  
    y_label="Diff=(CrystalFE-HomogeneousFE)/atom" 
    y_plot="DiffFE"
elif args.ycfev:  
    y_axis=thisdata[:,5]    
    y_label="Crystal Free Energy/volume"
    y_plot="cFEv"
elif args.yfv:  
    y_axis=thisdata[:,6]     
    y_label="Fraction of vacancies (fv)"
    y_plot="fv"
elif args.ygw:   
    y_axis=thisdata[:,7]
    y_label="Width of Gaussian (gwidth)"
    y_plot="gw"
elif args.ycol:   
    y_axis=thisdata[:,args.ycol]     
    y_label=args.ylab
    y_plot=args.ylab

plot_name=data_directory+"/plot_"+y_plot+"vs"+x_plot+"_"+fixed_quantity+fixed_value+".png"
plot_title=y_label+" vs "+x_label+" at Fixed "+fixed_quantity+"="+fixed_value

#Plot x-axis vs y-axis
plt.plot(x_axis, y_axis, color="purple")
plt.title(plot_title)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.savefig(plot_name)

plt.show()


