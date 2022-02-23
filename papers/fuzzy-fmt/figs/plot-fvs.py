#!/usr/bin/python3

#This program creates a plot of Free Energy difference vs gw at a specified 
#temperature and density from data in kT*n*alldat.dat (or kT*n*alldat_tensor.dat) 
#files which are generated as output data files by figs/new-melting.cpp
#The program compute-phasediagram-analysis.py is used to run new-melting 
#for various temperatures and densities.
#There are 2 plot options: 
#Plot option #1: FEdiff-vs-gw for a set of fv values with 5 seeds each (default)
#ie:  ./figs/plot-fvs.py --kT 0.5 --n 1.05 data/phase-diagram-test-fv-mcerr3
#Plot option #2: FEdiff-vs-gw for one fv value with 5 seeds
#ie:  ./figs/plot-fvs.py --kT 0.5 --n 1.07 --fv 0 data/phase-diagram-test-fv-mcerr3

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./figs/plot-fvs.py --kT [temp] --n [density] 
#                 --fv [OPTIONAL: enter 0, 1e-1, 1e-2, 1e-3, or 1e-4]  
#                 directory [data/phase-diagram-test-mcerr3] 


import os, glob
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a plot of FEdiff vs gw at a specified temperature, density, and one or more fraction of vacancies for five seeds.')

parser.add_argument('--kT', metavar='temperature', type=float, required=True,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', metavar='density', type=float, required=True,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--fv', metavar='fraction of vacancies', type=float,
                    help='fv - OPTIONAL enter 0, 1e-1, 1e-2, 1e-3, or 1e-4')
parser.add_argument('directory', metavar='directory', type=str,
                    help='directory with data to plot ie. data/phase-diagram-test-mcerr3 or data/phase-diagram-test-mcerr4') 
# parser.add_argument('--show', action='store_true',
                    # help='OPTIONAL show plot') 

args=parser.parse_args()

kT=args.kT
n=args.n
fv=args.fv

if fv is not None:
  doallfvs=0
else:
  doallfvs=1

fvs = (0, 1e-9, 1e-8, 1e-6, 1e-4, 1e-2)
seeds=(1,2,3,4,5)

#Plots all seeds for all fvs:
#if doallfvs == 1:
lowest_fe_difference =0
count_0=0
count_1=0
maxcount = len(fvs) * len(seeds)
maxfvcount = len(seeds) *2   #Two gw values
maxseedcount = len(seeds)
number_of_nans_at_fv = [maxfvcount,maxfvcount,maxfvcount,maxfvcount,maxfvcount,maxfvcount]
number_of_nans_at_fv_0 = [maxseedcount,maxseedcount]
number_of_nans_at_fv_1 = [maxseedcount,maxseedcount]
number_of_nans_at_fv_2 = [maxseedcount,maxseedcount]
number_of_nans_at_fv_3 = [maxseedcount,maxseedcount]
number_of_nans_at_fv_4 = [maxseedcount,maxseedcount]
number_of_nans_at_fv_5 = [maxseedcount,maxseedcount]
print (maxcount)

for fv in fvs:
  for seed in seeds:
    gw = []
    fe_difference = []  
    files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.9f*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed)))) 
    for i in range(len(files)):
      #print (i)
      #print(files[i])
      if i == 0:
         count_0=count_0+1  #only one gw has a file
      if i == 1:
         count_1=count_1+1  #both gw have files
    for f in files:
      data = np.loadtxt(f)
      gw.append(data[3])
      fe_difference.append(data[6])
      if fv == 1e-2:
        plt.plot(gw, fe_difference, '.-', color='yellow')
        number_of_nans_at_fv[5] = number_of_nans_at_fv[5]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_5[0] = number_of_nans_at_fv_5[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_5[1] = number_of_nans_at_fv_5[1]-1
      if fv == 1e-4:
        plt.plot(gw, fe_difference, '.-', color='blue')
        number_of_nans_at_fv[4] = number_of_nans_at_fv[4]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_4[0] = number_of_nans_at_fv_4[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_4[1] = number_of_nans_at_fv_4[1]-1
      if fv == 1e-6:
        plt.plot(gw, fe_difference, '.-', color='green')
        number_of_nans_at_fv[3] = number_of_nans_at_fv[3]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_3[0] = number_of_nans_at_fv_3[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_3[1] = number_of_nans_at_fv_3[1]-1
      if fv == 1e-8:
        plt.plot(gw, fe_difference, '.-', color='purple')
        number_of_nans_at_fv[2] = number_of_nans_at_fv[2]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_2[0] = number_of_nans_at_fv_2[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_2[1] = number_of_nans_at_fv_2[1]-1
      if fv == 1e-9:
        plt.plot(gw, fe_difference, '.-', color='orange')
        number_of_nans_at_fv[1] = number_of_nans_at_fv[1]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_1[0] = number_of_nans_at_fv_1[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_1[1] = number_of_nans_at_fv_1[1]-1
      if fv == 0:
        plt.plot(gw, fe_difference, '.-', color='red')
        number_of_nans_at_fv[0] = number_of_nans_at_fv[0]-1
        if data[3] == 0.01:
            number_of_nans_at_fv_0[0] = number_of_nans_at_fv_0[0]-1
        if data[3] == 0.05:
            number_of_nans_at_fv_0[1] = number_of_nans_at_fv_0[1]-1
      mcerror= data[12]	
      #print ('fe_difference=%g, gw=%g, fv=%g, seed=%g' %(data[6], data[3], fv, seed))
      if data[6] < lowest_fe_difference:
         lowest_fe_difference = data[6]
         lowest_fv = fv
         lowest_seed=seed 
         lowest_gw = data[3]     
plt.plot([],[], color='yellow', label='fv=1e-2')
plt.plot([],[], color='blue', label='fv=1e-4')
plt.plot([],[], color='green', label='fv=1e-6')
plt.plot([],[], color='purple', label='fv=1e-8')
plt.plot([],[], color='orange', label='fv=0.1e-9')
plt.plot([],[], color='red', label='fv=0')
plt.plot([],[], color='white', label='5 seeds for each fv')

plt.legend(loc='best')      
plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
plt.axhline(0, color='black')
plt.title('Free Energy difference vs gw  for kT=%g, n=%g, mcerror=%g' % (kT, n, mcerror))
plt.ylabel('FEdiff')
plt.xlabel('gw')
plt.show()

print(number_of_nans_at_fv)
number_of_nans = (maxcount - count_0)*2 + (maxcount-count_1)*1
print('count_0=%g, count_1=%g' % (count_0, count_1))
print('number_of_nans = %g' % (number_of_nans))
# print ('kT=%g, n=%g, lowest_gw=%g, lowest_fv=%g, lowest_seed=%g, lowest_fe_difference=%g' % (kT, n, lowest_gw, lowest_fv, lowest_seed, lowest_fe_difference))



plt.figure('Frequency of NANs')
plt.bar([1, 2, 3, 4, 5, 6], number_of_nans_at_fv)
plt.xticks([1, 2, 3, 4, 5, 6], ('0', '1e-9', '1e-8', '1e-6', '0.0001', '0.01'))
plt.title('Lowest gw=%g, at fv=%g and seed=%g' % (lowest_gw, lowest_fv, lowest_seed))
plt.ylabel('number of NANs')
plt.xlabel('fv')
plt.show()

#plt.figure('Frequency of NANs at various gw')
labels = ['0.01', '0.05']
x = np.arange(len(labels))
print(x)
width=0.1
fig, ax = plt.subplots()
rects1 = ax.bar(x-2.5*width, number_of_nans_at_fv_0, width, color='blue', label = 'fv=0')
rects2 = ax.bar(x-1.5*width, number_of_nans_at_fv_1, width, color='green', label = 'fv=1e-9')
rects3 = ax.bar(x-0.5*width, number_of_nans_at_fv_2, width, color='red', label = 'fv=1e-8')
rects4 = ax.bar(x+0.5*width, number_of_nans_at_fv_3, width, color='cyan', label = 'fv=1e-6')
rects5 = ax.bar(x+1.5*width, number_of_nans_at_fv_4, width, color='magenta', label = 'fv=0.0001')
rects6 = ax.bar(x+2.5*width, number_of_nans_at_fv_5, width, color='orange', label = 'fv=0.01')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_title('Lowest gw=%g, at fv=%g and seed=%g' % (lowest_gw, lowest_fv, lowest_seed))
ax.set_ylabel('number of NANs')
ax.set_xlabel('gw')

def autolabel(rects, colortext):
     for rect in rects:
        height = rect.get_height()		
        ax.annotate('{}'.format(height), xy=(rect.get_x() +rect.get_width()/2, height), 
        xytext=(0,3), textcoords = "offset points", ha = 'center', va='bottom', color=colortext)
autolabel(rects1, 'blue')
autolabel(rects2, 'green')
autolabel(rects3, 'red')
autolabel(rects4, 'cyan')
autolabel(rects5, 'magenta')
autolabel(rects6, 'orange')
ax.legend(loc='best') 
plt.show()

# plt.figure()
# #Plots all seeds at one fv per plot:
# #if doallfvs == 0:
# for fv in fvs:	
   # for seed in seeds:
      # gw = []
      # fe_difference = []  
      # files = sorted(list(glob.glob('%s/kT%.3f_n%.3f_fv%.6f*seed%g-alldat_tensor.dat' % (args.directory, kT, n, fv, seed)))) 
      # for f in files:
        # data = np.loadtxt(f)
        # gw.append(data[3])
        # fe_difference.append(data[6])
        # if seed == 1:
          # plt.plot(gw, fe_difference, '.-', color='blue')
        # if seed == 2:
          # plt.plot(gw, fe_difference, '.-', color='green')
        # if seed == 3:
          # plt.plot(gw, fe_difference, '.-', color='purple')
        # if seed == 4:
          # plt.plot(gw, fe_difference, '.-', color='orange')
        # if seed == 5:
          # plt.plot(gw, fe_difference, '.-', color='red')
        # mcerror= data[12]	
   # plt.plot([],[], color='white', label='fv=%g' % (fv))
   # plt.plot([],[], color='blue', label='seed1')
   # plt.plot([],[], color='green', label='seed2')
   # plt.plot([],[], color='purple', label='seed3')
   # plt.plot([],[], color='orange', label='seed4')
   # plt.plot([],[], color='red', label='seed5') 
   # plt.legend(loc='best')      
   # plt.axhspan(0.2, -0.2, color='black', alpha=0.15, lw=0)
   # plt.axhspan(0.02, -0.02, color='green', alpha=0.15, lw=0)
   # plt.axhline(0, color='black')
   # plt.title('Free Energy difference vs gw  for kT=%g, n=%g, mcerror=%g' % (kT, n, mcerror))
   # plt.ylabel('FEdiff')
   # plt.xlabel('gw')
   # plt.show()
   

