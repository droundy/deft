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



if len(sys.argv) < 2:
        print("Usage: %s TEMPERATURE" % sys.argv[0])
        exit(1)
temp=float(sys.argv[1])



################################## START INITIALIZATION #######################################
#                                                                                             #
#                                                                                             #
sigma = 2     #Sphere diameter
k_B = 1       #Boltzman's constant

numdata=30    #Number of numdensity data points
max_fillingfraction_handled = 0.55
sphere_volume = (sigma**3*np.pi/6)
numdensity = plt.linspace(0.0008,max_fillingfraction_handled/sphere_volume,numdata)   #number density

num_init = 1000
num_right = 200
num_left = 200
num_mid = 1000

datadir = 'data13'
try:
        os.mkdir(datadir)
except:
        pass
fload = datadir+'/fit_T%.5f_f'%temp
fsave = datadir+'/fit_T%.5f_f'%temp

print(temp)


def RG_first_pass(T,n,i):
        fnaught = RG.SWfid(T,n) + RG.SWfhs(T,n) + RG.a2(n)/k_B/T*n # SW (and Hughes) a2/kT is the same as Forte's f2
        f = fnaught
        return f + RG_dfi(n)

def RG_later_pass(n):
        f = f01_ext(n)
        return f + RG_dfi(n)

def RG_dfi(n):
        maxn = (max_fillingfraction_handled+0.0001)/sphere_volume
        # warning: Forte defines x as a density, we define it
        # as a dimensionless quantity that scales the density.
        maxx = np.minimum(1,maxn/n-1)

        if abs(maxx) < 1e-42:
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
        T=temp
        # eqn (5) from Forte 2011:
        IDvalue = ID2(n, maxx)
        IDvalueStar = ID2star(n, maxx)
        return -k_B*T*(IDvalue-IDvalueStar)/RG.VD(fn) # eqn (7), Forte 2011

def firstPass():
        global f01interp
        f01interp=lambda n: RG.fiterative(temp,n,0)
        data=[]
        t = time.time()
        lastprint = t
        for i in range(numdata):
            n = numdensity[i]
            print "%d of %d: "%(i,numdata),
            free_energy = RG_first_pass(temp,n,1)
            data.append([n,free_energy])
        elapsed = time.time() - t
        print(elapsed)

        np.savetxt(fsave,data)

def laterPass():
        f01_load()
        data=[]
        f02 = []
        for i in range(numdata):
                print "%d of %d: "%(i,numdata),
                n = numdensity[i]
                free_energy = RG_later_pass(n)
                f02.append(free_energy)
                data.append([n,free_energy])

        
        np.savetxt(fsave,data)


fn=1
while os.path.isfile(fload+'%d.out'%fn):
        fn+=1
fload+='%d.out'%(fn-1)
fsave+='%d.out'%fn
print('fn:',fn)


#                                                                                             #
#                                                                                             #
############################### END INITIALIATION #############################################








########################### START PLOTTING STUFF ##############################################        
#                                                                                             #
#                                                                                             #






def f0(n):
        return RG.fiterative(temp,n,0)


def f01_load():
        global f01interp
        f01data=np.loadtxt(fload)
        f01 = [f01data[i][1] for i in range(0,len(f01data))]  
        numdensity = [f01data[i][0] for i in range(0,len(f01data))]
        f01interp = interp1d(numdensity,f01,kind='cubic')



##integrand_xs = []
##integrand_args = []

def ID2(n, maxx):
##        This is $I_D$ from Forte's paper.
        

##        global integrandIDlistx
##        global integrandIDlistarg
        
        
        maxpower = -1e99
        dx = maxx/num_init
        xpoints = np.arange(dx/2,maxx, dx)
        kmax=0
        for k in range(len(xpoints)):
                if maxpower<onlyPower(n,xpoints[k],fn):
                        maxpower=onlyPower(n,xpoints[k],fn)
                        kmax=k
                
        #########FIND LEFT SIDE
        x_left = xpoints[kmax]
        k_left = kmax
        while k_left>0 and integrand_ID2(n,maxpower,x_left)>0.1:
                k_left -= 1
                x_left = xpoints[k_left]
        #########FIND RIGHT SIDE
        x_right = xpoints[kmax]
        k_right = kmax
        while k_right<len(xpoints)-1 and integrand_ID2(n,maxpower,x_right)>0.1:
                k_right += 1
                x_right = xpoints[k_right]
        

        integral_left = 0

        dx_left = x_left/num_left
        xpoints_left = np.arange(dx_left/2,x_left,dx_left)
        max_power_left = -1e99
        for k in range(len(xpoints_left)):
                max_power_left = max(max_power_left,onlyPower(n,xpoints_left[k],fn))
        

        integral_right = 0
 
        dx_right = (maxx-x_right)/num_right
        xpoints_right = np.arange(x_right+dx_right/2,maxx,dx_right)
        max_power_right = -1e99
        for k in range(len(xpoints_right)):
                max_power_right = max(max_power_right,onlyPower(n,xpoints_right[k],fn))


        
        integral_mid = 0
        
        dx_mid = (x_right-x_left)/num_mid
        xpoints_mid = np.arange(x_left+dx_mid/2,x_right,dx_mid)
        max_power_mid = -1e99
        for k in range(len(xpoints_mid)):
                max_power_mid = max(max_power_mid,onlyPower(n,xpoints_mid[k],fn))
        

        if max_power_left > max_power_mid:
                x_left = 0
                integral_left = 0                
        else:
                integral_left += integrate.midpoint(lambda x: integrand_ID2(n,max_power_left,x),0,x_left,num_left)*n
        if max_power_right > max_power_mid:
                x_right = maxx
                integral_right = 0
        else:
                integral_right += integrate.midpoint(lambda x: integrand_ID2(n,max_power_right,x),x_right,maxx,num_right)*n
        
        if x_left != 0 and x_right != maxx:
                integral_mid += integrate.midpoint(lambda x: integrand_ID2(n,max_power_mid,x),x_left,x_right,num_mid)*n

        else:
                dx_mid = (x_right-x_left)/num_mid
                xpoints_mid = np.arange(x_left+dx_mid/2,x_right,dx_mid)
                max_power_mid = -1e99
                for k in range(len(xpoints_mid)):
                        max_power_mid = max(max_power_mid,onlyPower(n,xpoints_mid[k],fn))
                integral_mid += integrate.midpoint(lambda x: integrand_ID2(n,max_power_mid,x),x_left,x_right,num_mid)*n
                

        if integral_left == 0 and integral_right == 0:
                print "TTm=%.2E   \r"%(max_power_mid),
                return np.log(integral_mid)+max_power_mid
        if integral_left != 0 and integral_right == 0:
                print "FTlm=%.2E   \r"%(max_power_left-max_power_mid),
                return np.log(integral_left*np.exp(max_power_left-max_power_mid)+integral_mid)+max_power_mid
        if integral_left == 0 and integral_right != 0:
                print "TFrm=%.2E   \r"%(max_power_right-max_power_mid),
                return np.log(integral_right*np.exp(max_power_right-max_power_mid)+integral_mid)+max_power_mid
        
                
        print "FFlm=%.2E FFrm=%.2E   \r"%(max_power_left-max_power_mid,max_power_right-max_power_mid),
        return np.log(integral_left*np.exp(max_power_left-max_power_mid)+integral_right*np.exp(max_power_right-max_power_mid)+integral_mid)+max_power_mid

        


def onlyPower(n,x,i):
        return (-RG.VD(i)/k_B/temp*(fbarD(temp,n,x,i) + ubarD(temp,n,x,i)))



##integrandIDlistn = []
##integrandIDlistx = []
##integrandIDlistarg = []
def integrand_ID2(n,maxpower,x):
        argument = np.exp(-maxpower-RG.VD(fn)/k_B/temp*(fbarD(temp,n,x,fn) + ubarD(temp,n,x,fn)))
##        integrandIDlistn.append(n)
##        integrandIDlistx.append(x)
##        integrandIDlistarg.append(argument)
        return argument



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
        if numdensity > 0.0008 and numdensity < max_fillingfraction_handled/sphere_volume:
                return f01interp(numdensity)
        return RG.fiterative(temp,numdensity,0)



        
        


        


def onlyPowerStar(n,x,i):
        return (-RG.VD(i)/k_B/temp*fbarD(temp,n,x,i))




        
        


def ID2star(n, maxx):

        maxpower = -1e99
        dx = maxx/num_init
        xpoints = np.arange(dx/2,maxx, dx)
        kmax=0
        for k in range(len(xpoints)):
                if maxpower<onlyPowerStar(n,xpoints[k],fn):
                        maxpower = onlyPowerStar(n,xpoints[k],fn)
                        kmax=k
                
        #########FIND LEFT SIDE
        x_left = xpoints[kmax]
        k_left = kmax
        while k_left>0 and integrand_ID2star(n,maxpower,x_left)>0.1:
                k_left -= 1
                x_left = xpoints[k_left]
        #########FIND RIGHT SIDE
        x_right = xpoints[kmax]
        k_right = kmax
        while k_right<len(xpoints)-1 and integrand_ID2star(n,maxpower,x_right)>0.1:
                k_right += 1
                x_right = xpoints[k_right]
        

        integral_left = 0

        dx_left = x_left/num_left
        xpoints_left = np.arange(dx_left/2,x_left,dx_left)
        max_power_left = -1e99
        for k in range(len(xpoints_left)):
                max_power_left = max(max_power_left,onlyPowerStar(n,xpoints_left[k],fn))
        

        integral_right = 0
 
        dx_right = (maxx-x_right)/num_right
        xpoints_right = np.arange(x_right+dx_right/2,maxx,dx_right)
        max_power_right = -1e99
        for k in range(len(xpoints_right)):
                max_power_right = max(max_power_right,onlyPowerStar(n,xpoints_right[k],fn))


        
        integral_mid = 0
        
        dx_mid = (x_right-x_left)/num_mid
        xpoints_mid = np.arange(x_left+dx_mid/2,x_right,dx_mid)
        max_power_mid = -1e99
        for k in range(len(xpoints_mid)):
                max_power_mid = max(max_power_mid,onlyPowerStar(n,xpoints_mid[k],fn))
        

        if max_power_left > max_power_mid:
                x_left = 0
                integral_left = 0                
        else:
                integral_left += integrate.midpoint(lambda x: integrand_ID2star(n,max_power_left,x),0,x_left,num_left)*n
        if max_power_right > max_power_mid:
                x_right = maxx
                integral_right = 0
        else:
                integral_right += integrate.midpoint(lambda x: integrand_ID2star(n,max_power_right,x),x_right,maxx,num_right)*n
        
        if x_left != 0 and x_right != maxx:
                integral_mid += integrate.midpoint(lambda x: integrand_ID2star(n,max_power_mid,x),x_left,x_right,num_mid)*n

        else:
                dx_mid = (x_right-x_left)/num_mid
                xpoints_mid = np.arange(x_left+dx_mid/2,x_right,dx_mid)
                max_power_mid = -1e99
                for k in range(len(xpoints_mid)):
                        max_power_mid = max(max_power_mid,onlyPowerStar(n,xpoints_mid[k],fn))
                integral_mid += integrate.midpoint(lambda x: integrand_ID2star(n,max_power_mid,x),x_left,x_right,num_mid)*n
                

        if integral_left == 0 and integral_right == 0:
                #print "TTm=%.2E   \r"%(max_power_mid),
                return np.log(integral_mid)+max_power_mid
        if integral_left != 0 and integral_right == 0:
                #print "FTlm=%.2E   \r"%(max_power_left-max_power_mid),
                return np.log(integral_left*np.exp(max_power_left-max_power_mid)+integral_mid)+max_power_mid
        if integral_left == 0 and integral_right != 0:
                #print "TFrm=%.2E   \r"%(max_power_right-max_power_mid),
                return np.log(integral_right*np.exp(max_power_right-max_power_mid)+integral_mid)+max_power_mid
        
                
        #print "FFlm=%.2E FFrm=%.2E   \r"%(max_power_left-max_power_mid,max_power_right-max_power_mid),
        return np.log(integral_left*np.exp(max_power_left-max_power_mid)+integral_right*np.exp(max_power_right-max_power_mid)+integral_mid)+max_power_mid

 
        
def integrand_ID2star(n,maxpower,x):        
        argument = np.exp(-maxpower-RG.VD(fn)/k_B/temp*fbarD(temp,n,x,fn))
        return argument        
        


if fn == 1:
        firstPass()
        print('First pass done')
        exit(0)

laterPass()



##plt.figure()
##plt.title('integrand vs x')
##plt.ylabel('integrand')
##plt.xlabel('x')
##for i in range(len(integrand_xs)):
##        plt.plot(integrand_xs[i], integrand_args[i])
##
##plt.show()



#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
