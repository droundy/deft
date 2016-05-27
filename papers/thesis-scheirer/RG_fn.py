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
numdata=100
numdensity = plt.linspace(0.0001,.2,numdata)

temp=float(sys.argv[1])
#temp=0.6
sigma=2
k_B=1

fload = 'data05/fit_T%.3f_f'%temp
fsave = 'data05/fit_T%.3f_f'%temp

print(temp)
f01interp=lambda n: 1


def RG_first_pass(T,n,i):
        # eqn (55) Forte 2011:
        global maxx
        global maxn
        maxn = 1/(sigma**3*np.pi/6)
        maxx = np.minimum(1,maxn/(n+1e-42)-1)
        fnaught = RG.SWfid(T,n) + RG.SWfhs(T,n) + RG.a2(n)/k_B/T*n # SW (and Hughes) a2/kT is the same as Forte's f2
        f = fnaught
        IDvalue = ID2(n)
        ID_refvalue = ID2star(n)
        dfi = -k_B*T*(IDvalue-ID_refvalue)/RG.VD(1) # eqn (7), Forte 2011
        f += dfi
        return f



def firstPass():
        global f01interp
        f01interp=lambda n: RG.fiterative(temp,n,0)
        data=[]
        y = plt.linspace(0,1,len(numdensity))
        t = time.time()
        for i in range(0,len(y)):
                print "%d of %d     \r"%(i,len(y)),
                y[i]=RG_first_pass(temp,numdensity[i],1)
                data.append([numdensity[i],y[i]])
        elapsed = time.time() - t
        print(elapsed)        

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
                maxx = np.minimum(1,maxn/(n+1e-42)-1)
                
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
        maxpower = -1e99
        numdensityMark10=(plt.linspace(0.00001,maxx*1.05,100))
        xtemp=numdensityMark10[0]
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPower(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPower(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]
        numdensityMark10=plt.linspace(xtemp*0.9,xtemp*1.1,100)
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPower(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPower(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]
        
        numdensityMark10=plt.linspace(xtemp*0.98,xtemp*1.02,1000)
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPower(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPower(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]
        
        #maxpower = sp.optimize.fmin(lambda x: -onlyPowerExperimental(n,x,maxx,fn),maxpower)

        #integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),maxx*0.1,maxx*0.9)[0]*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx*0.1,100)*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),maxx*0.9,maxx,100)*n

        #integral=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx,1000)*n
        integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),0,maxx)[0]*n
        return np.log(integral)+maxpower  #returns log of ID
        #return integral


def onlyPower(n,x,i):
        return (-RG.VD(i)/k_B/temp*(fbarD(temp,n,x,i) + ubarD(temp,n,x,i)))

def onlyPowerExperimental(n,x,maxx,fn):
        if(x>=maxx): return -1e99
        if x<=0 : return -1e99
        return onlyPower(n,x,fn)

#integrandIDlistn = []
integrandIDlistx = []
integrandIDlistarg = []
def integrand_ID2(n,maxpower,x):        
        argument = np.exp(-maxpower-RG.VD(fn)/k_B/temp*(fbarD(temp,n,x,fn) + ubarD(temp,n,x,fn)))
        #integrandIDlistn.append(n)
        integrandIDlistx.append(x)
        integrandIDlistarg.append(argument)
        return argument




def fbarD(T,n,x,i):
        iplusx = f01_ext(n*(1+x))
        iminusx = f01_ext(n*(1-x))
        nochangex = f01_ext(n)
        value = (iplusx + iminusx)/2 - nochangex
        return value

# Averge scaled potential
# eqn (11), Forte 2011
def ubarD(T,n,x,i):
        value = (RG.u(temp,n*(1+x),0,i) + RG.u(temp,n*(1-x),0,i))/2 - RG.u(temp,n,0,i)
        return value











        

def f01_ext(numdensity):
        if numdensity > 0.0001 and numdensity < 0.2:
                return f01interp(numdensity)
        return RG.fiterative(temp,numdensity,0)



        
        


        


def onlyPowerStar(n,x,i):
        return (-RG.VD(i)/k_B/temp*fbarD(temp,n,x,i))

def onlyPowerStarExperimental(n,x,maxx,fn):
        if(x>=maxx): return -1e99
        if x<=0 : return -1e99
        return onlyPowerStar(n,x,fn)



        
        


def ID2star(n):

        maxpower = -1e99
        xtemp=numdensity[0]
        numdensityMark10=(plt.linspace(0.00001,maxx*1.05,100))
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPowerStar(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPowerStar(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]
        
        numdensityMark10=plt.linspace(xtemp*0.9,xtemp*1.1,100)
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPowerStar(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPowerStar(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]

        numdensityMark10=plt.linspace(xtemp*0.98,xtemp*1.02,1000)
        for k in range(0,len(numdensityMark10)):
                if numdensityMark10[k] > maxx: break
                #maxpower = max(maxpower,onlyPower(n,numdensityMark10[k],fn))
                if onlyPowerStar(n,numdensityMark10[k],fn)>maxpower:
                        maxpower=onlyPowerStar(n,numdensityMark10[k],fn)
                        xtemp=numdensityMark10[k]
        #maxpower2 = sp.optimize.fmin(lambda x: -onlyPowerStar(n,x,fn),0)
        #maxpower = sp.optimize.fmin(lambda x: -onlyPowerStarExperimental(n,x,maxx,fn),maxpower)
        integral = sp.integrate.quad(lambda x: integrand_ID2star(n,maxpower,x),0,maxx)[0]*n
        #print('scipy integral:',integral)
        #integral=integrate.midpoint(lambda x: integrand_ID2star(n,maxpower,x),0,maxx,1000)*n

        #integral = sp.integrate.quad(lambda x: integrand_ID2star(n,maxpower,x),maxx*0.1,maxx*0.9)[0]*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2star(n,maxpower,x),0,maxx*0.1,100)*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2star(n,maxpower,x),maxx*0.9,maxx,100)*n

        return np.log(integral)+maxpower  #returns log of ID
        #return integral


        
        '''maxpower = -1e10
        for k in range(0,len(numdensity)):
                if numdensity[k] > maxx: break
                maxpower = max(maxpower,onlyPowerStar(n,numdensity[k],fn))
        integral=integrate.midpoint(lambda x: integrand_ID2star(n,maxpower,x),0,maxx,1000)*n
        return np.log(integral)+maxpower  #returns log of ID
        #return integral'''


        
        
def integrand_ID2star(n,maxpower,x):        
        argument = np.exp(-maxpower-RG.VD(fn)/k_B/temp*fbarD(temp,n,x,fn))
        return argument        
        






def fit2(n):
        T=temp
        f = f01_ext(n)
        # eqn (5) from Forte 2011:
        IDvalue = ID2(n)
        IDvalueStar = ID2star(n)
        dfi = -k_B*T*(IDvalue-IDvalueStar)/RG.VD(fn) # eqn (7), Forte 2011
        f += dfi

        return f










#################################################
#       CALL THE PLOTS YOU WANT TO PLOT         #
#                                               #
#testload()
#f01_load()
#testf01()
if fn == 1:
        firstPass()
        sys.exit("First Pass Done")

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

#np.savetxt('integrandIDlistn',integrandIDlistn)
#np.savetxt('integrandIDlistx',integrandIDlistx)
#np.savetxt('integrandIDlistarg',integrandIDlistarg)

#plt.figure()
#plt.title('integrand vs n')
#plt.ylabel('integrand')
#plt.xlabel('number density')
#plt.xlim([0.015,0.030])
#plt.ylim([-0.01,1e112])
#plt.plot(integrandIDlistx,integrandIDlistarg)
#plt.savefig('meeting/23may2016/integrand_int100_03_T%.3f.png'%temp)
#plt.show()

plt.figure()
plt.title('integrand vs x')
plt.ylabel('integrand')
plt.xlabel('x')
plt.plot(integrandIDlistx,integrandIDlistarg)
plt.show()



#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
