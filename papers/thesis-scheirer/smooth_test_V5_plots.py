from __future__ import division
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import numpy as np
import matplotlib
import pylab as plt
import os
import sys
import RG
import SW
import time


temp = plt.linspace(0.6,1.28,20)

data1 = np.loadtxt('RG_cotangent_data1.out')
data2 = np.loadtxt('RG_cotangent_data2.out')
data3 = np.loadtxt('RG_cotangent_data3.out')
data4 = np.loadtxt('RG_cotangent_data4.out')
data5 = np.loadtxt('RG_cotangent_data5.out')
data6 = np.loadtxt('RG_cotangent_data6.out')
datasnaft = np.loadtxt('snaft2.out')



x = plt.linspace(1e-20,.2,4000)                         # x-axis grid space (used for when plotting number density on x-axis)
xmin=(3/(4.0*np.pi))*(1e-20)*(1/(SW.R)**3)
xmax=(3/(4.0*np.pi))*(0.2)*(1/(SW.R)**3)
xff = plt.linspace(xmin,xmax,4000)


nLs = [datasnaft[i][1] for i in range(0,len(datasnaft))]           # list of the left number density solutions obtained from fsolve
nRs = [datasnaft[i][2] for i in range(0,len(datasnaft))]           # list of the right number density solutions obtained from fsolve
Tlists = [datasnaft[i][0] for i in range(0,len(datasnaft))]        # list of the corresponding temperatures for which the above number densitites were found

ffLs = [i*((4*np.pi*(SW.R)**3)/3) for i in nLs]           # converts left number density to filling fraction
ffRs = [i*((4*np.pi*(SW.R)**3)/3) for i in nRs]           # converts right number density to filling fraction


nL1 = [data1[i][1] for i in range(0,len(data1))]           # list of the left number density solutions obtained from fsolve
nR1 = [data1[i][2] for i in range(0,len(data1))]           # list of the right number density solutions obtained from fsolve
Tlist1 = [data1[i][0] for i in range(0,len(data1))]        # list of the corresponding temperatures for which the above number densitites were found

ffL1 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL1]           # converts left number density to filling fraction
ffR1 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR1]           # converts right number density to filling fraction


nL2 = [data2[i][1] for i in range(0,len(data2))]           # list of the left number density solutions obtained from fsolve
nR2 = [data2[i][2] for i in range(0,len(data2))]           # list of the right number density solutions obtained from fsolve
Tlist2 = [data2[i][0] for i in range(0,len(data2))]        # list of the corresponding temperatures for which the above number densitites were found

ffL2 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL2]           # converts left number density to filling fraction
ffR2 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR2]           # converts right number density to filling fraction


nL3 = [data3[i][1] for i in range(0,len(data3))]           # list of the left number density solutions obtained from fsolve
nR3 = [data3[i][2] for i in range(0,len(data3))]           # list of the right number density solutions obtained from fsolve
Tlist3 = [data3[i][0] for i in range(0,len(data3))]        # list of the corresponding temperatures for which the above number densitites were found

ffL3 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL3]           # converts left number density to filling fraction
ffR3 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR3]           # converts right number density to filling fraction


nL4 = [data4[i][1] for i in range(0,len(data4))]           # list of the left number density solutions obtained from fsolve
nR4 = [data4[i][2] for i in range(0,len(data4))]           # list of the right number density solutions obtained from fsolve
Tlist4 = [data4[i][0] for i in range(0,len(data4))]        # list of the corresponding temperatures for which the above number densitites were found

ffL4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL4]           # converts left number density to filling fraction
ffR4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR4]           # converts right number density to filling fraction


nL5 = [data5[i][1] for i in range(0,len(data5))]           # list of the left number density solutions obtained from fsolve
nR5 = [data5[i][2] for i in range(0,len(data5))]           # list of the right number density solutions obtained from fsolve
Tlist5 = [data5[i][0] for i in range(0,len(data5))]        # list of the corresponding temperatures for which the above number densitites were found

ffL5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL5]           # converts left number density to filling fraction
ffR5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR5]           # converts right number density to filling fraction


nL6 = [data6[i][1] for i in range(0,len(data6))]           # list of the left number density solutions obtained from fsolve
nR6 = [data6[i][2] for i in range(0,len(data6))]           # list of the right number density solutions obtained from fsolve
Tlist6 = [data6[i][0] for i in range(0,len(data6))]        # list of the corresponding temperatures for which the above number densitites were found

ffL6 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL6]           # converts left number density to filling fraction
ffR6 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR6]           # converts right number density to filling fraction



def liq_vap_Tvsff():
  plt.figure()
  #plt.plot(ffL,Tlist,color='#f36118',linewidth=2)
  #plt.plot(ffR,Tlist,color='c',linewidth=2)

  plt.plot(ffLs,Tlists,color='c',linewidth=2)
  plt.plot(ffRs,Tlists,color='c',linewidth=2)
  
  plt.plot(ffL1,Tlist1,color='r',linewidth=2)
  plt.plot(ffR1,Tlist1,color='r',linewidth=2)

  plt.plot(ffL2,Tlist2,color='g',linewidth=2)
  plt.plot(ffR2,Tlist2,color='g',linewidth=2)

  plt.plot(ffL3,Tlist3,color='b',linewidth=2)
  plt.plot(ffR3,Tlist3,color='b',linewidth=2)

  plt.plot(ffL4,Tlist4,color='k',linewidth=2)
  plt.plot(ffR4,Tlist4,color='k',linewidth=2)

  #plt.plot(ffL5,Tlist5,color='c',linewidth=2)
  #plt.plot(ffR5,Tlist5,color='c',linewidth=2)

  #plt.plot(ffL6,Tlist6,color='r',linewidth=2)
  #plt.plot(ffR6,Tlist6,color='r',linewidth=2)

  
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  plt.xlim(-.05,max(ffR1)+.05)
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
  #plt.show()
  plt.savefig('meeting/13may2016/Tvsff_RG.png'%temp)



liq_vap_Tvsff()






