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
test_dn = (max_fillingfraction_handled/(4*np.pi/3))/numdata
#numdensity = np.arange(0.0001,max_fillingfraction_handled/(4*np.pi/3),test_dn)
numdensity = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),numdata)

if len(sys.argv) < 2:
        print("Usage: %s TEMPERATURE" % sys.argv[0])
        exit(1)
temp=float(sys.argv[1])

sigma=2
k_B=1

datadir = 'data09'
try:
        os.mkdir(datadir)
except:
        pass
fload = datadir+'/fit_T%.3f_f'%temp
fsave = datadir+'/fit_T%.3f_f'%temp

print(temp)
#f01interp=lambda n: 1


def RG_first_pass(T,n,i):
        # eqn (55) Forte 2011:
        global maxx
        global maxn
        maxn = (max_fillingfraction_handled+0.01)/(sigma**3*np.pi/6)
        '''NEED TO ADD SMALL THING ABOVE IF FAIL'''
        maxx = np.minimum(1,maxn/n-1)
        
        if abs(maxx) < 1e-42:
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
                print 'maxx is zero'
                
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
        lastprint = t
        for i in range(0,len(y)):
                print "%d of %d     \r"%(i,len(y)),
                '''if time.time() - lastprint > 5:
                        print "%d of %d     \n"%(i,len(y)),
                        lastprint = time.time()'''
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
        #numdensity2 = np.arange(0.0001,max_fillingfraction_handled/(4*np.pi/3),test_dn)
        numdensity2 = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),numdata)
        #numdensity2 = plt.linspace(0.0001,.2,numdata)
        f01_load()
        data=[]
        f02=[]
        for i in range(0,len(numdensity2)):
                n = numdensity2[i]
                maxn = (max_fillingfraction_handled+0.01)/(sigma**3*np.pi/6)
                # warning: Forte defines x as a density, we define it
                # as a dimensionless quantity that scales the density.
                maxx = np.minimum(1,maxn/n-1)

                if abs(maxx) < 1e-42:
                        print 'maxx is zero'
                        print 'maxx is zero'
                        print 'maxx is zero'
                        print 'maxx is zero'
                        print 'maxx is zero'

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


'''integrand_xs = []
integrand_args = []'''

def ID2(n):
        ''' This is $I_D$ from Forte's paper.
        '''

        '''global integrandIDlistx
        global integrandIDlistarg'''
        
        
        maxpower = -1e99
        dx = maxx/1000
        #maxxx = maxx+dx
        #xpoints = plt.linspace(dx/2.0,maxx,1000)
        xpoints = np.arange(dx/2,maxx, dx)
        kmax=0
        for k in range(len(xpoints)):
                if maxpower<onlyPower(n,xpoints[k],fn):
                        maxpower=onlyPower(n,xpoints[k],fn)
                        kmax=k
                #maxpower = max(maxpower, onlyPower(n,xpoints[k],fn))

        xLEFT = xpoints[kmax]
        kLeft=kmax
        while kLeft>0 and integrand_ID2(n,maxpower,xLEFT)>0.1:
                kLeft-=1
                xLEFT=xpoints[kLeft]
        xRIGHT= xpoints[kmax]
        kRight=kmax
        while kRight<len(xpoints)-1 and integrand_ID2(n,maxpower,xRIGHT)>0.1:
                kRight+=1
                xRIGHT=xpoints[kRight]
        '''integrandIDlistx = []
        integrandIDlistarg = []'''

        #integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),maxx*0.1,maxx*0.9)[0]*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx*0.1,100)*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),maxx*0.9,maxx,100)*n

        integralLEFT = 0

        dxLeft=xLEFT/200
        xpointsLeft=np.arange(dxLeft/2,xLEFT,dxLeft)
        maxPowerLeft=-1e99
        for k in range(len(xpointsLeft)):
                maxPowerLeft=max(maxPowerLeft,onlyPower(n,xpointsLeft[k],fn))
        integralLEFT += integrate.midpoint(lambda x: integrand_ID2(n,maxPowerLeft,x),0,xLEFT,200)*n

        '''integrand_xs.append(integrandIDlistx)
        integrand_args.append(integrandIDlistarg)
        integrandIDlistx = []
        integrandIDlistarg = []'''




        integralRIGHT = 0

        if xRIGHT > maxx:
                xRIGHT = maxx
                maxPowerRight=0
        else:
                dxRight=(maxx-xRIGHT)/200
                xpointsRight=np.arange(xRIGHT+dxRight/2,maxx,dxRight)
                maxPowerRight=-1e99
                for k in range(len(xpointsRight)):
                        maxPowerRight=max(maxPowerRight,onlyPower(n,xpointsRight[k],fn))
                integralRIGHT += integrate.midpoint(lambda x: integrand_ID2(n,maxPowerRight,x),xRIGHT,maxx,200)*n

        '''integrand_xs.append(integrandIDlistx)
        integrand_args.append(integrandIDlistarg)
        integrandIDlistx = []
        integrandIDlistarg = []'''


        integralMID = 0
        
        dxMid=(xRIGHT-xLEFT)/1000
        xpointsMid=np.arange(xLEFT+dxMid/2,xRIGHT,dxMid)
        maxPowerMid=-1e99
        for k in range(len(xpointsMid)):
                maxPowerMid=max(maxPowerMid,onlyPower(n,xpointsMid[k],fn))
        integralMID += integrate.midpoint(lambda x: integrand_ID2(n,maxPowerMid,x),xLEFT,xRIGHT,1000)*n


        
        
                

        #integral=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx,1000)*n
        #print "maxpower:", maxpower

        #print 'maxpower:', maxPowerLeft+maxPowerRight+maxPowerMid

        
        #integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),0,maxx)[0]*n
        
        '''integrand_xs.append(integrandIDlistx)
        integrand_args.append(integrandIDlistarg)
        integrandIDlistx = []
        integrandIDlistarg = []'''
        return np.log(integralLEFT*np.exp(maxPowerLeft-maxPowerMid)+integralRIGHT*np.exp(maxPowerRight-maxPowerMid)+integralMID)+maxPowerMid

        #return np.log(integral)+maxPowerLeft+maxPowerRight+maxPowerMid  #returns log of ID
        #return integral


def onlyPower(n,x,i):
        return (-RG.VD(i)/k_B/temp*(fbarD(temp,n,x,i) + ubarD(temp,n,x,i)))

def onlyPowerExperimental(n,x,maxx,fn):
        if(x>=maxx): return -1e99
        if x<=0 : return -1e99
        return onlyPower(n,x,fn)

#integrandIDlistn = []
'''integrandIDlistx = []
integrandIDlistarg = []'''
def integrand_ID2(n,maxpower,x):
        argument = np.exp(-maxpower-RG.VD(fn)/k_B/temp*(fbarD(temp,n,x,fn) + ubarD(temp,n,x,fn)))
        #integrandIDlistn.append(n)
        '''integrandIDlistx.append(x)
        integrandIDlistarg.append(argument)'''
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
        if numdensity > 0.0001 and numdensity < max_fillingfraction_handled/(4*np.pi/3):
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
        dx = maxx/1000
        #maxxx = maxx+dx
        #xpoints = plt.linspace(dx/2.0,maxx,1000)
        xpoints = np.arange(dx/2,maxx, dx)
        kmax=0
        for k in range(len(xpoints)):
                if maxpower<onlyPowerStar(n,xpoints[k],fn):
                        maxpower=onlyPowerStar(n,xpoints[k],fn)
                        kmax=k
                #maxpower = max(maxpower, onlyPower(n,xpoints[k],fn))

        xLEFT = xpoints[kmax]
        kLeft=kmax
        while kLeft>0 and integrand_ID2star(n,maxpower,xLEFT)>0.1:
                kLeft-=1
                xLEFT=xpoints[kLeft]
        xRIGHT= xpoints[kmax]
        kRight=kmax
        while kRight<len(xpoints)-1 and integrand_ID2star(n,maxpower,xRIGHT)>0.1:
                kRight+=1
                xRIGHT=xpoints[kRight]

        #integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),maxx*0.1,maxx*0.9)[0]*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx*0.1,100)*n
        #integral+=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),maxx*0.9,maxx,100)*n

        integralLEFT = 0

        dxLeft=xLEFT/200
        xpointsLeft=np.arange(dxLeft/2,xLEFT,dxLeft)
        maxPowerLeft=-1e99
        for k in range(len(xpointsLeft)):
                maxPowerLeft=max(maxPowerLeft,onlyPowerStar(n,xpointsLeft[k],fn))
        integralLEFT += integrate.midpoint(lambda x: integrand_ID2star(n,maxPowerLeft,x),0,xLEFT,200)*n






        integralRIGHT = 0

        if xRIGHT > maxx:
                xRIGHT = maxx
                maxPowerRight=0
        else:
                dxRight=(maxx-xRIGHT)/200
                xpointsRight=np.arange(xRIGHT+dxRight/2,maxx,dxRight)
                maxPowerRight=-1e99
                for k in range(len(xpointsRight)):
                        maxPowerRight=max(maxPowerRight,onlyPowerStar(n,xpointsRight[k],fn))
                integralRIGHT += integrate.midpoint(lambda x: integrand_ID2star(n,maxPowerRight,x),xRIGHT,maxx,200)*n


        integralMID = 0

        dxMid=(xRIGHT-xLEFT)/1000
        xpointsMid=np.arange(xLEFT+dxMid/2,xRIGHT,dxMid)
        maxPowerMid=-1e99
        for k in range(len(xpointsMid)):
                maxPowerMid=max(maxPowerMid,onlyPowerStar(n,xpointsMid[k],fn))
        integralMID += integrate.midpoint(lambda x: integrand_ID2star(n,maxPowerMid,x),xLEFT,xRIGHT,1000)*n



        
                

        #integral=integrate.midpoint(lambda x: integrand_ID2(n,maxpower,x),0,maxx,1000)*n
        #print "maxpower:", maxpower
        #print 'maxpower:', maxPowerLeft+maxPowerRight+maxPowerMid
        #integral = sp.integrate.quad(lambda x: integrand_ID2(n,maxpower,x),0,maxx)[0]*n
        return np.log(integralLEFT*np.exp(maxPowerLeft-maxPowerMid)+integralRIGHT*np.exp(maxPowerRight-maxPowerMid)+integralMID)+maxPowerMid
        #return integral

        
        
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

'''plt.figure()
plt.title('integrand vs x')
plt.ylabel('integrand')
plt.xlabel('x')
for i in range(len(integrand_xs)):
        plt.plot(integrand_xs[i], integrand_args[i])

plt.show()'''



#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
