import scipy as sp
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import pylab as plt
import matplotlib
import RG
import SW
import numpy as np
import time
import integrate
import sys
import os


###############################################################################################
# Author: Ryan Scheirer                                                                       #
# Email: scheirer@oregonstate.edu                                                             #
# Date: February 2016                                                                         #
                                                                                              #
# Uses fsolve to find the common tangent of the free energy density vs number density...      #
# ...this then constructs the temp vs filling fraction liquid-vapor coexistence plot, total...#
# ...grand free energy per volume, and many more fun plots.                                   #
###############################################################################################





################################## START INITIALIZATION #######################################
#                                                                                             #
#                                                                                             #
### Normal temperature linspace (useful for quick troubleshooting) ###
#temp = plt.linspace(0.4,1.3,40)
numdensity = plt.linspace(0.0001,.2,1000)
numdensity2 = plt.linspace(0.0001,.2,1000)
temp = plt.linspace(.6,1.28,20)
temp=temp[0]
#temp=float(sys.argv[1])
sigma=2
k_B=1
#                                                                                             #
#                                                                                             #
############################### END INITIALIATION #############################################









########################### START PLOTTING STUFF ##############################################        
#                                                                                             #
#                                                                                             #





def f01_load():
        global f01interp
        fname = 'data/fit_T%.3f_f26.out'%temp
        f01data=np.loadtxt(fname)
        f01 = [f01data[i][1] for i in range(0,len(f01data))]  
        numdensity = [f01data[i][0] for i in range(0,len(f01data))]
        f01interp = interp1d(numdensity,f01,kind='cubic')

        #numdensity_interp=np.linspace(0.0001,0.2,100)
        #plt.plot(numdensity,f01,'o',numdensity_interp,f01interp(numdensity_interp))
        #plt.show()



def f02_load():
        global f02interp
        fname = 'data/fit_T%.3f_f28.out'%temp
        f02data=np.loadtxt(fname)
        f02 = [f02data[i][1] for i in range(0,len(f02data))]  
        numdensity = [f02data[i][0] for i in range(0,len(f02data))]
        f02interp = interp1d(numdensity,f02,kind='cubic')

        #numdensity_interp=np.linspace(0.0001,0.2,100)
        #plt.plot(numdensity,f01,'o',numdensity_interp,f01interp(numdensity_interp))
        #plt.show()


def f03_load():
        global f03interp
        fname = 'data/fit_T%.3f_f30.out'%temp
        f03data=np.loadtxt(fname)
        f03 = [f03data[i][1] for i in range(0,len(f03data))]  
        numdensity = [f03data[i][0] for i in range(0,len(f03data))]
        f03interp = interp1d(numdensity,f03,kind='cubic')

def f04_load():
        global f04interp
        fname = 'data/fit_T%.3f_f32.out'%temp
        f04data=np.loadtxt(fname)
        f04 = [f04data[i][1] for i in range(0,len(f04data))]  
        numdensity = [f04data[i][0] for i in range(0,len(f04data))]
        f04interp = interp1d(numdensity,f04,kind='cubic')

def f05_load():
        global f05interp
        fname = 'data/fit_T%.3f_f33.out'%temp
        f05data=np.loadtxt(fname)
        f05 = [f05data[i][1] for i in range(0,len(f05data))]  
        numdensity = [f05data[i][0] for i in range(0,len(f05data))]
        f05interp = interp1d(numdensity,f05,kind='cubic')




def f00_tot(n):
        return RG.ftot(temp,n,0)

def f01_tot(n):
        return f01_ext(n) + RG.a1SW(n)*n


def f02_tot(n):
        return f02_ext(n) + RG.a1SW(n)*n

def f03_tot(n):
        return f03_ext(n) + RG.a1SW(n)*n
def f04_tot(n):
        return f04_ext(n) + RG.a1SW(n)*n
def f05_tot(n):
        return f05_ext(n) + RG.a1SW(n)*n




def plotstuff():
        numdensity3 = plt.linspace(0.0001,.2,1000)
        f01_load()
        f02_load()
        f03_load()
        f04_load()
        f05_load()
        ysw=[]
        yrg0=[]
        yrg1=[]
        yrg2=[]
        yrg3=[]
        yrg4=[]
        yrg5=[]
        print (f02interp(0.1))
        for i in range(0,len(numdensity3)):
                ysw.append(SW.ftot(temp,numdensity3[i]))
                yrg0.append(f00_tot(numdensity3[i]))
                yrg1.append(f01_tot(numdensity3[i])-ysw[i])
                yrg2.append(f02_tot(numdensity3[i])-ysw[i])
                yrg3.append(f03_tot(numdensity3[i])-ysw[i])
                yrg4.append(f04_tot(numdensity3[i])-ysw[i])
                yrg5.append(f05_tot(numdensity3[i])-ysw[i])
                


        plt.figure()
        plt.title('free energy difference')
        plt.ylabel('free energy difference')
        plt.xlabel('number density')
        #plt.plot(numdensity,ysw,color='r',linewidth=6)
        #plt.plot(numdensity,yrg0,'b-',linewidth=3)
        plt.plot(numdensity3,yrg1,'g--',linewidth=2)
        plt.plot(numdensity3,yrg2,'r--',linewidth=2)
        plt.plot(numdensity3,yrg3,color='k',linewidth=2)
        plt.plot(numdensity3,yrg4,color='b',linewidth=2)
        plt.plot(numdensity3,yrg5,color='g',linewidth=2)
        #plt.show()
        
        plt.savefig('figs_diff/RG_difference_26_33_T%.3f.png'%temp)
        plt.close()


        yrg0=[]
        yrg1=[]
        yrg2=[]
        yrg3=[]
        yrg4=[]
        yrg5=[]

        for i in range(0,len(numdensity3)):
                #ysw.append(SW.ftot(temp,numdensity3[i]))
                yrg0.append(f00_tot(numdensity3[i]))
                yrg1.append(f01_tot(numdensity3[i]))
                yrg2.append(f02_tot(numdensity3[i]))
                yrg3.append(f03_tot(numdensity3[i]))
                yrg4.append(f04_tot(numdensity3[i]))
                yrg5.append(f05_tot(numdensity3[i]))

        plt.figure()
        plt.title('free energy difference')
        plt.ylabel('free energy difference')
        plt.xlabel('number density')
        #plt.plot(numdensity,ysw,color='r',linewidth=6)
        #plt.plot(numdensity3,yrg0,'r',linewidth=2)
        #plt.plot(numdensity3,yrg1,'g--',linewidth=2)
        #plt.plot(numdensity3,yrg2,'r--',linewidth=2)
        #plt.plot(numdensity3,yrg3,color='b',linewidth=2)
        #plt.plot(numdensity3,yrg4,color='g',linewidth=2)
        plt.plot(numdensity3,yrg5,color='g',linewidth=1)
        #plt.show()
        plt.savefig('RG_FREE.pdf')
        #plt.close()



        

def f01_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f01interp(numdensity)
        return RG.fiterative(temp,numdensity,0)


def f02_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f02interp(numdensity)
        return RG.fiterative(temp,numdensity,0)

def f03_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f03interp(numdensity)
        return RG.fiterative(temp,numdensity,0)

def f04_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f04interp(numdensity)
        return RG.fiterative(temp,numdensity,0)

def f05_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f05interp(numdensity)
        return RG.fiterative(temp,numdensity,0)

        
        

plotstuff()
 



#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
