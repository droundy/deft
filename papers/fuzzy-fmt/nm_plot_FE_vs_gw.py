#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./nm_plot_FE_vs_gw.py [filename.dat]
#Look for a filename.dat (file with data to plot) with a name like:  
#FE_vs_gw_kT2_n1.3_fv0_dx0.5_mcerror0.001_mcconstant5_mcprefactor50000_seeds10.dat


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

data_length = len(thisdata)

#show constants on plot:
dx=thisdata[0,0]
mc_error=thisdata[0,1]
mc_constant=thisdata[0,10]
mc_prefactor=thisdata[0,11]
fv=thisdata[0,13]
kT=thisdata[0,14]
n=thisdata[0,15]

#variables:
FE = thisdata[:,7]
gw = thisdata[:,12]

#gw_old=gw[0]
#FE_sum=FE[0]
#j=1
#mean_FE={}

#for i in range(1,data_length):
   #if (gw[i] == gw_old):
      #FE_sum=FE_sum+FE[i]
      #j=j+1
   #else :
      ##print FE_sum
      ##print j
      #mean_FE[gw[i-1]]=FE_sum/j
      ##print mean_FE
      #gw_old=gw[i] 
      #j=1
      #FE_sum=FE[i]
      ##print
##print FE_sum
##print j
#mean_FE[gw[i]]=FE_sum/j 
##print mean_FE

gw_old=gw[0]
FE_allvalues_at_gw=[1]
mean_FE_uncertainty={}
mean_FE={}
FE_allvalues_at_gw[0]=thisdata[0,7]
j=1

for i in range(1,data_length):
   if (gw[i] == gw_old):
      FE_allvalues_at_gw.append(thisdata[i,7])
      j=j+1
   else :
      #print FE_allvalues_at_gw  #debug
      #print j  #debug
      mean_FE_uncertainty[gw[i-1]]=np.std(FE_allvalues_at_gw)/np.sqrt(j)
      mean_FE[gw[i-1]]=np.mean(FE_allvalues_at_gw)
      #print mean_FE_uncertainty  #debug
      gw_old=gw[i] 
      FE_allvalues_at_gw=[1] #re-initialize
      FE_allvalues_at_gw[0]=thisdata[i,7]
      j=1
#print FE_allvalues_at_gw  #debug
#print j  #debug
mean_FE_uncertainty[gw[i]]=np.std(FE_allvalues_at_gw)/np.sqrt(j-1)
#print mean_FE_uncertainty  #debug
mean_FE[gw[i]]=np.mean(FE_allvalues_at_gw)
#print mean_FE  #debug
#print  #debug
#print  #debug

#print 'TEST', mean_FE
#mean_FE_bruteforce=(thisdata[0,7]+thisdata[1,7]+thisdata[2,7]+thisdata[3,7]+thisdata[4,7]+thisdata[5,7]+thisdata[6,7]+thisdata[7,7]+thisdata[8,7]+thisdata[9,7])/10
#print 'Compare to brute force', mean_FE_bruteforce

for gw in mean_FE.keys():
   print '#gw     mean_FE[gw]     mean_FE_uncertainty[gw]   kT   n  fv'  
   print gw, mean_FE[gw], mean_FE_uncertainty[gw], kT, n, fv #for data file
   #TO DO: create data out file with the values of gw, mean_FE[gw], mean_FE_uncertainty[gw]
   plt.errorbar(gw, mean_FE[gw], mean_FE_uncertainty[gw], color='blue', fmt='o')
   #The errorbars are small - ZOOM in to see! 
   
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')

plt.title('FE vs gw at kT %g n %g fv %g' % (kT, n, fv))
plt.ylabel('FE')
plt.xlabel('gw (with dx=%g mc_error=%g mc_constant=%g mc_prefactor=%g)' % (dx, mc_error, mc_constant, mc_prefactor))

plt.savefig(args.filedat+".png" )
plt.show()


