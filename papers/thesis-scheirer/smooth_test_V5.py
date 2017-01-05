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



### Initialize stuff ###
numdensity = plt.linspace(0.0001,.2,1000)
numdensity2 = plt.linspace(0.0001,.2,1000)
T = plt.linspace(.6,1.28,20)

Lguess = .0015  # Initial left "n" guess 
Rguess = 0.1    # Initial right "n" guess

data = []

for p in range(2,20):

        temp=T[p]




        def fload():
                # Loads RG data and splits the file into two lists for number density and RG free energy
                # Also uses the interp1d function from scipy to interpolate between data points
                global finterp
                fname = 'data/fit_T%.3f_f6.out'%temp
                fdata=np.loadtxt(fname)
                f = [fdata[i][1] for i in range(0,len(fdata))]  
                numdensity = [fdata[i][0] for i in range(0,len(fdata))]
                finterp = interp1d(numdensity,f,kind='cubic')



        def f_tot(n):
                #  Returns total free energy data with the missing a1SW(n) term that was not used in RG process
                return f_ext(n) + RG.a1SW(n)*n


        def f_ext(numdensity):
                # Fsolve might need to go outside of our minimum or maximum density that RG was used with, thus this function...
                # ... sets the values of the free energy outside those bounds.
                if numdensity > 0.0001 and numdensity < 0.2:
                        return finterp(numdensity)
                return RG.fiterative(temp,numdensity,0)






        def frgp_ext(numdensity):
                # Fsolve might need to go outside of our minimum or maximum density that RG was used with, thus this function...
                # ... sets the values of the free energy outside those bounds. 
                if numdensity > 0.0001 and numdensity < 0.2:
                        return frgpinterp(numdensity)
                if numdensity > 0.2:
                        return frgpinterp(0.2)
                return frgpinterp(0.0001)
        

        def dfrgp_dn(n):
            ### Derivative of total free energy per volume w.r.t. n  (NOTE: at fixed volume and pressure, this is mu) ###
            dn = 1e-6*n
            return (frgp_ext(n + dn) - frgp_ext(n - dn))/(2*dn)










        def plotstuff():
                global rgp
                global frgpinterp
                numdensity3 = plt.linspace(0.0001,.2,1000)
                fload()
                ysw=[]  # List for SnAFT SW fluid total free energy data
                yrgdiffsw=[]  # List for fitted free energy difference (total RG free energy - total SnAFT SW free energy) 
                savdiff=[]  # List for difference between raw free energy difference and fitted data from savgo
                rgfittot=[]  # List for fitted RG total free energy
                rg=[]  # List for RG total free energy
                rgdiff = []  # List for difference between RG total free energy and the savgo fitted total free energy
                for i in range(0,len(numdensity3)):
                        ysw.append(SW.ftot(temp,numdensity3[i]))
                        yrgdiffsw.append(f_tot(numdensity3[i])-ysw[i]) 
 
    

                savgo = savgol_filter(yrgdiffsw,1,0)  #  Uses Savitzky Golay filter from scipy package to fit data (gets rid of noise)
                # Use the arguments of (function,1,0) if you wish to not apply the filter

                for j in range(0,len(numdensity3)):
                        # Builds the fitted RG total free energy (rgfittot)                
                        # Builds difference between raw free energy difference and fitted data from savgo
                        # Not necessiary, I was just curious how much the savgo filter adjusted the data
                        rgfittot.append(SW.ftot(temp,numdensity3[j])+savgo[j])
                        savdiff.append(yrgdiffsw[j]-savgo[j])
                        rg.append(f_tot(numdensity3[j]))
                        rgdiff.append(rgfittot[j]-rg[j])

                
                frgpinterp = interp1d(numdensity,rgfittot,kind='cubic')
              
                plt.figure()
                plt.title('savgo - rg')
                plt.ylabel('free energy')
                plt.xlabel('number density')
                #plt.xlim([0,0.11])
                #plt.ylim([-0.01,0.02])


        
                #plt.plot(numdensity3,savgo,color='red',linewidth=2)
                #plt.plot(numdensity3,sav,color='red',linewidth=2)
        
        
                #plt.plot(numdensity3,yrgdiffsw,color='green',linewidth=2)
                #plt.plot(numdensity3,savgo,color='red',linewidth=1)
                plt.plot(numdensity3,rgdiff,color='green',linewidth=1)
                yrgp = []
                dyrgp = []
                madyrgp = []
                for i in range(0,len(numdensity3)):
                        yrgp.append(frgp_ext(numdensity3[i]))
                        dyrgp.append(dfrgp_dn(numdensity3[i]))
                        if numdensity3[i] < 0.085 and numdensity3[i]>0.001:
                                madyrgp.append(dyrgp[i])
                madyrgp = sum(madyrgp)/len(madyrgp)
                #print madyrgp
                cotangentGuess = []
                for i in range(0,len(numdensity3)):
                        cotangentGuess.append(yrgp[i]-madyrgp*numdensity3[i])
                #plt.plot(numdensity3,yrgp)
                #plt.plot(numdensity3,rg,color='r',linewidth=2)
                #plt.plot(numdensity3,dyrgp)
                #plt.plot(numdensity3,cotangentGuess)
                #plt.plot(0.00079008,frgp_ext(0.00079008),'ko')
                #plt.plot(0.09041596,frgp_ext(0.09041596),'ro')
                #plt.plot(0.00262281,frgp_ext(0.00262281)-madyrgp*0.00262281,'ko')
                #plt.plot(0.08894228,frgp_ext(0.08894228)-madyrgp*0.08894228,'ro')
                #plt.plot(0.00079008,frgp_ext(0.00079008)-madyrgp*0.00079008,'ko')
                #plt.plot(0.09041596,frgp_ext(0.09041596)-madyrgp*0.09041596,'ro')
                #plt.plot(numdensity3,ysw,color='b',linewidth=2)
                #plt.plot(numdensity3,rgfittot, color='green',linewidth=2)
                #plt.show()
                #plt.savefig('meeting/13may2016/rgdiff_5_w1_po0_T%.3f.png'%temp)
                



        plotstuff()





        def g(x):
                ### Just a general function which includes conditions for finding the cotangent
                x1=x[0]                                 # left guess "n"
                x2=x[1]                                 # right guess "n"
                y1=frgp_ext(x[0])                      # f(x1)
                y2=frgp_ext(x[1])                      # f(x2)
                dy1=dfrgp_dn(x[0])            # df_dx(x1)
                dy2=dfrgp_dn(x[1])            # df_dx(x2)

                             
                out=[(dy1-dy2)]                         # Condition 1: df_dx(x1) = df_dx(x2)
                out.append(dy1-((y2-y1)/(x2-x1)))       # Condition 2: df_dx(x1) = slope between the two positions
                return out



        sol = fsolve(g,[Lguess,Rguess])         # Magic happens here
        Lguess = sol[0]
        Rguess = sol[1]                 # Resets the next right guess number density for fsolve as the previous right solution
        data.append([temp,sol[0],sol[1]])
        print(sol)
        print(p)
        #print(data)
        np.savetxt('RG_cotangent_data6.out',data)












