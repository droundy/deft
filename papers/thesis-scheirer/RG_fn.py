from __future__ import division

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
numdata=1000
max_fillingfraction_handled = 0.84
numdensity = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),numdata)

if len(sys.argv) < 2:
        print("Usage: %s TEMPERATURE" % sys.argv[0])
        exit(1)
temp=float(sys.argv[1])
sigma=2
k_B=1

fload = 'data/fit_T%.3f_f'%temp
fsave = 'data/fit_T%.3f_f'%temp



def firstPass():
        data=[]
        y = plt.linspace(0,1,len(numdensity))
        t = time.time()
        for i in range(0,len(y)):
                print "%d of %d     \r"%(i,len(y)),
                y[i]=RG.fiterative(temp,numdensity[i],1)
                data.append([numdensity[i],y[i]])
        elapsed = time.time() - t
        print(elapsed)        

        np.savetxt(fsave,data)


fn=1
while os.path.isfile(fload+'%d.out'%fn):
        fn+=1
fload+='%d.out'%(fn-1)
fsave+='%d.out'%fn
if fn == 1:
        firstPass()
        print("First Pass Done")
        


#                                                                                             #
#                                                                                             #
############################### END INITIALIATION #############################################








########################### START PLOTTING STUFF ##############################################        
#                                                                                             #
#                                                                                             #




def testfit02():
        global f02
        global numdensity2
        global maxn
        global maxx
        numdensity2 = plt.linspace(0.0001,.2,numdata)
        f01_load()
        data=[]
        f02=[]
        for i in range(0,len(numdensity2)):
                n = numdensity2[i]
                maxn = 1/(sigma**3*np.pi/6)
                # warning: Forte defines x as a density, we define it
                # as a dimensionless quantity that scales the density.
                maxx = np.minimum(1.0, maxn/(n+1e-42)-1)
                print "%d of %d     \r"%(i,len(numdensity2)),
                r = float(numdensity2[i])
                f02.append(float(fit2(r)))
                data.append([numdensity2[i],f02[i]])
                

        #plt.figure()
        #plt.plot(numdensity2,f02)
        #plt.show()
        #print(f01tot[0:20])

        
        np.savetxt(fsave,data)


def f0(n):
        return RG.fiterative(temp,n,0)


def f01_load():
        global f01interp
        f01data=np.loadtxt(fload)
        f01 = [f01data[i][1] for i in range(0,len(f01data))]  
        numdensity = [f01data[i][0] for i in range(0,len(f01data))]
        f01interp = interp1d(numdensity,f01,kind='cubic')

        #numdensity_interp=np.linspace(0.0001,0.2,100)
        #plt.plot(numdensity,f01,'o',numdensity_interp,f01interp(numdensity_interp))
        #plt.show()



def ID2(n):
        ''' This is $I_D$ from Forte's paper.
        '''
        maxpower = -1e10
        for k in range(0,len(numdensity)):
                if numdensity[k] > maxx: break
                maxpower = max(maxpower,onlyPower(n,numdensity[k],2))
        integral=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx,1000)*n
        return np.log(integral)+maxpower  #returns log of ID
        #return integral


def onlyPower(n,x,i):
        return (-RG.VD(i)/k_B/temp*(fbarD(temp,n,x,i) + ubarD(temp,n,x,i)))



def integrand_ID2(n,maxpower,x):
        return np.exp(-maxpower-RG.VD(2)/k_B/temp*(fbarD(temp,n,x,2) + ubarD(temp,n,x,2)))




def fbarD(T,n,x,i):
        iplusx = f01_ext(n*(1+x))
        iminusx = f01_ext(n*(1-x))
        nochangex = f01_ext(n)
        return (iplusx + iminusx)/2 - nochangex

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,n,x,i):
        return (RG.u(temp,n*(1+x),0,i) + RG.u(temp,n*(1-x),0,i))/2 - RG.u(temp,n,0,i)











        

def f01_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f01interp(numdensity)
        return RG.fiterative(temp,numdensity,0)



        
        


        


def onlyPowerStar(n,x,i):
        return (-RG.VD(i)/k_B/temp*fbarD(temp,n,x,i))




        
        


def ID2star(n):
        maxpower = -1e10
        for k in range(0,len(numdensity)):
                if numdensity[k] > maxx: break
                maxpower = max(maxpower,onlyPowerStar(n,numdensity[k],2))
        integral=integrate.midpoint(lambda x: integrand_ID2star(n,maxpower,x),0,maxx,1000)*n
        return np.log(integral)+maxpower  #returns log of ID
        #return integral


        
        
def integrand_ID2star(n,maxpower,x):        
        argument = np.exp(-maxpower-RG.VD(2)/k_B/temp*fbarD(temp,n,x,2))
        return argument        
        






def fit2(n):
        T=temp
        f = f01_ext(n)
        # eqn (5) from Forte 2011:
        IDvalue = ID2(n)
        IDvalueStar = ID2star(n)
        dfi = -k_B*T*(IDvalue-IDvalueStar)/RG.VD(2) # eqn (7), Forte 2011
        f += dfi

        return f










#################################################
#       CALL THE PLOTS YOU WANT TO PLOT         #
#                                               #
#testload()
#f01_load()
#testf01()
testfit02()
#april16th()
#cotangent_t0()
#expTest()
#cotangent_tloop()
#liq_vap_co_tvsff()
#cotangent_fsolve_breakdown()
#liq_vap_co_pvsT()
#                                               #
#                                               #
#################################################

#np.savetxt('figs/snaft2.out',data)

#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
