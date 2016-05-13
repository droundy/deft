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
import os
import sys


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
temp = plt.linspace(0.6,1.28,20)
#                                                                                             #
#                                                                                             #
############################### END INITIALIATION #############################################









########################### START PLOTTING STUFF ##############################################        
#                                                                                             #
#                                                                                             #




fns = []
fnames = []

def fns_load():
        global fnames
        dirname = os.getcwd()
        files = os.listdir(dirname+'/data')
        for arg in files:
                if 'f5.out' in arg:
                        fnames.append(dirname+'/data/'+arg)
                        

        for arg in fnames:
                f05data = np.loadtxt(arg)
                f05 = [f05data[i][1] for i in range(0,len(f05data))]
                numdensity = [f05data[i][0] for i in range(0,len(f05data))]
                f05interp = interp1d(numdensity,f05,kind='cubic')
                fns.append(f05interp)


def fns_tot(n,i):
        return fns_ext(n,i) + RG.a1SW(n)*n


def fns_ext(numdensity,i):
        if numdensity > 0.0001 and numdensity < 0.2:
                return fns[i](numdensity)
        return RG.fiterative(temp[i],numdensity,0)
        



def plotstuff():
        numdensity3 = plt.linspace(0.0001,.2,1000)
        fns_load()
        
        for i in range(0,len(fns)):
                y = []
                for j in range(0,len(numdensity3)):
                        y.append(float(fns_tot(float(numdensity3[j]),i)))

                plt.figure()
                plt.plot(numdensity3,y)
                plt.title(fnames[i])
                savename = os.getcwd()
                #savename += '/figs/'
                savename += fnames[i].split('data/')[1].split('.out')[0]
                savename += '.pdf'
                #print savename
                plt.savefig(savename)
                #plt.show()

        







plotstuff()
 





















#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
