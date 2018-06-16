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


#thisdata_n0_01_8000000 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC8000000.dat")
thisdata_n0_01_800000 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC800000.dat")
thisdata_n0_01_80000 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC80000.dat")
thisdata_n0_01_8000 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC8000.dat")
thisdata_n0_01_800 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC800.dat")
thisdata_n0_01_400 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC400.dat")
thisdata_n0_01_80 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC80.dat")
thisdata_n0_01_8 = np.loadtxt("DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC8.dat")


MC=[8]

#for i in range(1):
#for i in range(len(MC)):
#    print MC[i]
#    data_file="DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC"+str(MC[i])+".dat"
#    data=np.loadtxt(data_file)
#    print data
#    new_data=[MC[i], data[0,1]]
#    print new_data
#    new_array=[new_data]
#print new_array

new_array = []

MC=[8, 80, 400, 800, 8000, 80000, 800000]

for i in range(len(MC)):
#    print MC[i]
    data_file="DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC"+str(MC[i])+".dat"
    data=np.loadtxt(data_file)
#    print data
    new_data=[MC[i],data[0,1], data[0,2], data[0,3], data[0,4]]
#    print new_data
    new_array.append(new_data)
#print new_array



#if args.xdx:  
#    x_axis=thisdata[:,0]    
#    x_label="dx"
#    x_plot="dx"
#elif args.xcphi_1:  
#    x_axis=thisdata[:,1]     
#    x_label="cphi_1"
#    x_plot="c_phi_1"
#elif args.xcphi_2:  
#    x_axis=thisdata[:,2]     
#    x_label="cphi_2"
#    x_plot="cphi_2"
#elif args.xcphi3:   
#    x_axis=thisdata[:,3]     
#    x_label="cphi_3"
#    x_plot="cphi3"
#elif args.xcFE_per_vol:  
#    x_axis=thisdata[:,4]     
#    x_label="Crystal Free Energy per volume" 
#    x_plot="cFE_per_vol"
#elif args.xcol:   
#    x_axis=thisdata[:,args.xcol]     
#    x_label=args.xlab
#    x_plot=args.xlab
    
#if args.ydx:  
#    y_axis=thisdata[:,0]    
#    y_label="dx"
#    y_plot="dx"
#elif args.ycphi_1:  
#    y_axis=thisdata[:,1]     
#    y_label="cphi_1"
#    y_plot="c_phi_1"
#elif args.ycphi_2:  
#    y_axis=thisdata[:,2]     
#    y_label="cphi_2"
#    y_plot="cphi_2"
#elif args.ycphi3:   
#    y_axis=thisdata[:,3]     
#    y_label="cphi_3"
#    y_plot="cphi3"
#elif args.ycFE_per_vol:  
#    y_axis=thisdata[:,4]     
#    y_label="Crystal Free Energy per volume" 
#    y_plot="cFE_per_vol"
#elif args.ycol:   
#    y_axis=thisdata[:,args.ycol]     
#    y_label=args.ylab
#    y_plot=args.ylab

#plot_name=data_directory+"/plot_"+y_plot+"vs"+x_plot+"_"+fixed_quantity+fixed_value+".png"  
#plot_name="DATA_LOGS/PLOTS_n1.0/plot_"+y_plot+"vs"+x_plot+".png"   
#if args.ptname:
#    plot_name=data_directory+"/plot_"+y_plot+"vs"+x_plot+"_"+fixed_quantity+fixed_value+"_"+args.ptname+".png"


####plot_title=y_label+" vs "+x_label+" at Fixed "+fixed_quantity+"="+fixed_value
#plot_title=y_label+" vs "+x_label

#a=0
#plt.plot(x_axis, y_axis, color="purple")
#b=plt.gca().get_ylim()[0]-10

#Plot x-axis vs y-axis
#if b < 0 :
#    plt.axhspan(a, b, color='b', alpha=0.15, lw=0)
#else : plt.axhspan(a, -10, color='b', alpha=0.15, lw=0)
#plt.plot(x_axis, y_axis, color="purple")



cFEpervol=1  #plots cFEperVol vs dx for all MCpoints

cphi_1=1     #plots cphi_1 vs dx for all MCpoints
cphi_2=1
cphi_3=1

MC8=1     #plots cphi_1, cphi2 and cphi3 vs dx for specified MCpoints
MC80=1
MC400=1
MC800=1
MC8000=1
MC80000=1
MC800000=1
MC8000000=0

custom=0


if custom > 0 :
  MC=[800, 8000, 80000]   #CHOSE enteries
  for i in range(len(MC)):
    print MC[i]
    data_file="DATA_LOGS/PLOTS_n0_01/DATA_LOG_PLOT_n0_01_MC"+str(MC[i])+".dat"
    print data_file
    thisdata=np.loadtxt(data_file)
    plt.plot(thisdata[:,0], thisdata[:,1], label='MC'+str(MC[i]))   #CHOOSE COLUMES

  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_"+name_me    #FILL IN NAME

#  xlabel="dx"  #col 0			#SET xlabel and ylabel
#  ylabel="cphi1"  #col 1
#  ylabel="cphi2"  #col 2
#  ylabel="cphi3"  #col 3
#  ylabel="cFEperVol"  #col 4
  xlabel=""		
  ylabel=""

  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  plt.legend()
  plt.savefig(plot_name)


if cFEpervol > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,4], label = 'cFEperVol_MC8')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,4], label = 'cFEperVol_MC80')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,4], label = 'cFEperVol_MC400')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,4], label = 'cFEperVol_MC800')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,4], label = 'cFEperVol_MC8000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,4], label = 'cFEperVol_MC80000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,4], label = 'cFEperVol_MC800000')
  ##plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,4], label = 'cFEperVol_MC8000000') #won't plot one point?
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_cCEpervol"
  plot_title="cFEperVol vs dx"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cFEperVol")
  plt.legend()
  plt.savefig(plot_name)

if cphi_1 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,1], label = 'cphi_1_MC8')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,1], label = 'cphi_1_MC80')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,1], label = 'cphi_1_MC400')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,1], label = 'cphi_1_MC800')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,1], label = 'cphi_1_MC8000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,1], label = 'cphi_1_MC80000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,1], label = 'cphi_1_MC800000')
  ##plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,1], label = 'cphi_1_MC8000000') #won't plot one point?
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_cphi1"

  plot_title="cphi_1 vs dx"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi_1")	
  plt.legend()
  plt.savefig(plot_name)

if cphi_2 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,2], label = 'cphi_2_MC8')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,2], label = 'cphi_2_MC80')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,2], label = 'cphi_2_MC400')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,2], label = 'cphi_2_MC800')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,2], label = 'cphi_2_MC8000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,2], label = 'cphi_2_MC80000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,2], label = 'cphi_2_MC800000')
  ##plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,2], label = 'cphi_2_MC8000000') #won't plot one point?
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_cphi2"
  plot_title="cphi2 vs dx"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi_2")
  plt.legend()
  plt.savefig(plot_name)

if cphi_3 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,3], label = 'cphi_3_MC8')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,3], label = 'cphi_3_MC80')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,3], label = 'cphi_3_MC400')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,3], label = 'cphi_3_MC800')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,3], label = 'cphi_3_MC8000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,3], label = 'cphi_3_MC80000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,3], label = 'cphi_3_MC800000')
  ##plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,3], label = 'cphi_3_MC8000000') #won't plot one point?
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_cphi3"
  plot_title="cphi_3 vs dx"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi_3")
  plt.legend()
  plt.savefig(plot_name)

if MC8 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,1], label = 'cphi_1_MC8')
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,2], label = 'cphi_2_MC8')
  plt.plot(thisdata_n0_01_8[:,0], thisdata_n0_01_8[:,3], label = 'cphi_3_MC8')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC8_cphi123"
  plot_title="cphi123 vs dx with MC8"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC80 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,1], label = 'cphi_1_MC80')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,2], label = 'cphi_2_MC80')
  plt.plot(thisdata_n0_01_80[:,0], thisdata_n0_01_80[:,3], label = 'cphi_3_MC80')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC80_cphi123"
  plot_title="cphi123 vs dx with MC80"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC400 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,1], label = 'cphi_1_MC400')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,2], label = 'cphi_2_MC400')
  plt.plot(thisdata_n0_01_400[:,0], thisdata_n0_01_400[:,3], label = 'cphi_3_MC400')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC400_cphi123"
  plot_title="cphi vs dx with MC400"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC800 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,1], label = 'cphi_1_MC800')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,2], label = 'cphi_2_MC800')
  plt.plot(thisdata_n0_01_800[:,0], thisdata_n0_01_800[:,3], label = 'cphi_3_MC800')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC800_cphi123"
  plot_title="cphi123 vs dx with MC800"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC8000 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,1], label = 'cphi_1_MC8000')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,2], label = 'cphi_2_MC8000')
  plt.plot(thisdata_n0_01_8000[:,0], thisdata_n0_01_8000[:,3], label = 'cphi_3_MC8000')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC8000_cphi123"
  plot_title="cphi123 vs dx with MC8000"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC80000 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,1], label = 'cphi_1_MC80000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,2], label = 'cphi_2_MC80000')
  plt.plot(thisdata_n0_01_80000[:,0], thisdata_n0_01_80000[:,3], label = 'cphi_3_MC80000')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC80000_cphi123"
  plot_title="cphi123 vs dx with MC80000"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC800000 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,1], label = 'cphi_1_MC800000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,2], label = 'cphi_2_MC800000')
  plt.plot(thisdata_n0_01_800000[:,0], thisdata_n0_01_800000[:,3], label = 'cphi_3_MC800000')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC800000_cphi123"
  plot_title="cphi123 vs dx with MC800000"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

if MC8000000 > 0 :
  plt.figure()
  plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,1], label = 'cphi_1_MC8000000')
  plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,2], label = 'cphi_2_MC8000000')
  plt.plot(thisdata_n0_01_8000000[:,0], thisdata_n0_01_8000000[:,3], label = 'cphi_3_MC8000000')
  plot_name="DATA_LOGS/PLOTS_n0_01/plot_DATA_LOG_n0_01_MC8000000_cphi123"
  plot_title="cphi123 vs dx with MC 8000000"
  plt.title(plot_title)
  plt.xlabel("dx")
  plt.ylabel("cphi")
  plt.legend()
  plt.savefig(plot_name)

#plt.scatter(thisdata_n1_8[:,0], thisdata_n1_8[:,1], label = 'cphi_1_MC8')
#plt.scatter(thisdata_n1_80[:,0], thisdata_n1_80[:,1], label = 'cphi_1_MC80')
#plt.scatter(thisdata_n1_400[:,0], thisdata_n1_400[:,1], label = 'cphi_1_MC400')
#plt.scatter(thisdata_n1_800[:,0], thisdata_n1_800[:,1], label = 'cphi_1_MC800')
#plt.scatter(thisdata_n1_8000[:,0], thisdata_n1_8000[:,1], label = 'cphi_1_MC8000')
#plt.scatter(thisdata_n1_80000[:,0], thisdata_n1_80000[:,1], label = 'cphi_1_MC80000')
#plt.scatter(thisdata_n1_800000[:,0], thisdata_n1_800000[:,1], label = 'cphi_1_MC800000')
#plt.scatter(thisdata_n1_8000000[:,0], thisdata_n1_8000000[:,1], label = 'cphi_1_MC8000000')

#plt.scatter(thisdata_n1_800[:,0], thisdata_n1_800[:,2], label = 'cphi_2_MC800', color='blue')
#plt.scatter(thisdata_n1_800[:,0], thisdata_n1_800[:,3], label = 'cphi_3_MC800')

####plt.scatter(x_axis, y_axis)
#plt.title(plot_title)
#plt.xlabel("dx")
#plt.ylabel("cphi")
####plt.xlabel(x_label)
####plt.ylabel(y_label)
#plt.legend()


#plt.savefig(plot_name)

plt.show()


