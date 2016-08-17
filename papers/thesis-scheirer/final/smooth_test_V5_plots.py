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


#temp = plt.linspace(0.6,1.28,20)
#temp = np.concatenate((plt.linspace(1.20,1.21,6),plt.linspace(1.21+0.01/6,1.22,15),plt.linspace(1.22+0.01/15,1.23,6)),axis=0)


##numdata2=250    #Number of numdensity data points
##max_fillingfraction_handled = 0.15
##ffs = plt.linspace(0.0001,max_fillingfraction_handled,numdata2)
##temp = []
##
##for ff in ffs:
##	temp.append(1.215 - ((0.5)/(0.15)**3)*abs(ff-0.15)**3)
##


data1 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data1.out')
data2 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data2.out')
data3 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data3.out')
##data4 = np.loadtxt('RG_final_data/08_27Ts/cotangent_data/RG_cotangent_data4.out')
data5 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data5.out')
##data6 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data6.out')
##data7 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data7.out')
##data8 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data8.out')

datasnaft = np.loadtxt('snaft2.out')


numdata=1000
max_fillingfraction_handled = 0.55

x = plt.linspace(1e-20,max_fillingfraction_handled/(4*np.pi/3),numdata)                         # x-axis grid space (used for when plotting number density on x-axis)
xmin=(3/(4.0*np.pi))*(1e-20)*(1/(SW.R)**3)
xmax=(3/(4.0*np.pi))*(0.2)*(1/(SW.R)**3)
xff = plt.linspace(xmin,xmax,numdata)



  
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


##nL4 = [data4[i][1] for i in range(0,len(data4))]           # list of the left number density solutions obtained from fsolve
##nR4 = [data4[i][2] for i in range(0,len(data4))]           # list of the right number density solutions obtained from fsolve
##Tlist4 = [data4[i][0] for i in range(0,len(data4))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL4]           # converts left number density to filling fraction
##ffR4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR4]           # converts right number density to filling fraction
##
##
nL5 = [data5[i][1] for i in range(0,len(data5))]           # list of the left number density solutions obtained from fsolve
nR5 = [data5[i][2] for i in range(0,len(data5))]           # list of the right number density solutions obtained from fsolve
Tlist5 = [data5[i][0] for i in range(0,len(data5))]        # list of the corresponding temperatures for which the above number densitites were found

ffL5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL5]           # converts left number density to filling fraction
ffR5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR5]           # converts right number density to filling fraction


##nL6 = [data6[i][1] for i in range(0,len(data6))]           # list of the left number density solutions obtained from fsolve
##nR6 = [data6[i][2] for i in range(0,len(data6))]           # list of the right number density solutions obtained from fsolve
##Tlist6 = [data6[i][0] for i in range(0,len(data6))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL6 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL6]           # converts left number density to filling fraction
##ffR6 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR6]           # converts right number density to filling fraction
##
##
##nL7 = [data7[i][1] for i in range(0,len(data7))]           # list of the left number density solutions obtained from fsolve
##nR7 = [data7[i][2] for i in range(0,len(data7))]           # list of the right number density solutions obtained from fsolve
##Tlist7 = [data7[i][0] for i in range(0,len(data7))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL7 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL7]           # converts left number density to filling fraction
##ffR7 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR7]           # converts right number density to filling fraction
##
##
##
##nL8 = [data8[i][1] for i in range(0,len(data8))]           # list of the left number density solutions obtained from fsolve
##nR8 = [data8[i][2] for i in range(0,len(data8))]           # list of the right number density solutions obtained from fsolve
##Tlist8 = [data8[i][0] for i in range(0,len(data8))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL8 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL8]           # converts left number density to filling fraction
##ffR8 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR8]           # converts right number density to filling fraction
##



white_MD = np.loadtxt('white_MD.dat') # MD

forte_data = np.loadtxt('forte_data.dat') # not converged

forte_data_conv = np.loadtxt('forte_data_conv.dat') # Converged

mrhoWhite = (0.7-0.0)/(345-73)
brhoWhite = 0-mrhoWhite*73
mTWhite = (1.2-0.8)/(76-400)
bTWhite = 1.2-mTWhite*76

rhostar_MD = white_MD[:,0]*mrhoWhite + brhoWhite
T_MD = white_MD[:,1]*mTWhite + bTWhite
eta_MD = rhostar_MD*np.pi/6



mrho = (0.8-0.0)/(504-94)
brho = 0-mrho*94

mT = (0.65-1.85)/(345-52)
bT = 0.65-mT*345
rhostar_forte = forte_data[:,0]*mrho + brho
T_forte = forte_data[:,1]*mT + bT
eta_forte = rhostar_forte*np.pi/6
rhostar_conv = forte_data_conv[:,0]*mrho + brho
T_forte_conv = forte_data_conv[:,1]*mT + bT
eta_conv = rhostar_conv*np.pi/6



def liq_vap_Tvsff():
  plt.figure()
  #plt.plot(ffL,Tlist,color='#f36118',linewidth=2)
  #plt.plot(ffR,Tlist,color='c',linewidth=2)

  #plt.plot(eta_MD,T_MD,'yo',label='White 2000')
  
  plt.plot(eta_forte,T_forte,'b--',label='Forte 2011, non-converged RGT')


  plt.plot(eta_conv,T_forte_conv,'r--',label='Forte 2011, converged RGT')


  

  plt.plot(ffLs,Tlists,'c',linewidth=2)
  plt.plot(ffRs,Tlists,'c',linewidth=2)
  
##  plt.plot(ffL1,Tlist1,'ro',linewidth=2)
##  plt.plot(ffR1,Tlist1,'ro',linewidth=2)
##
##  plt.plot(ffL2,Tlist2,'go',linewidth=2)
##  plt.plot(ffR2,Tlist2,'go',linewidth=2)
##
##  plt.plot(ffL3,Tlist3,'bo',linewidth=2)
##  plt.plot(ffR3,Tlist3,'bo',linewidth=2)
##
##  plt.plot(ffL4,Tlist4,'ko',linewidth=2)
##  plt.plot(ffR4,Tlist4,'ko',linewidth=2)
##
  plt.plot(ffL5,Tlist5,'yo',linewidth=2)
  plt.plot(ffR5,Tlist5,'yo',linewidth=2)

##  plt.plot(ffL6,Tlist6,color='c',linewidth=2)
##  plt.plot(ffR6,Tlist6,color='c',linewidth=2)
##
##  plt.plot(ffL7,Tlist7,color='orange',linewidth=2)
##  plt.plot(ffR7,Tlist7,color='orange',linewidth=2)
##
##  plt.plot(ffL8,Tlist8,'co',linewidth=2)
##  plt.plot(ffR8,Tlist8,'co',linewidth=2)

  
  plt.xlabel(r'$filling fraction$')
  plt.ylabel(r'$temperature$')
  #plt.xlim(-.05,max(ffR1)+.05)
  #plt.xlim([0,0.40])
  #plt.ylim([0.8,1.39])
  plt.title('liquid-vapor coexistence '+r'$\lambda_{SW}=1.5$')
  plt.show()
  #plt.savefig('meeting/13may2016/Tvsff_RG_white2.png'%temp)



liq_vap_Tvsff()






