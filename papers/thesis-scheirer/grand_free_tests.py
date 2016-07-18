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



data1 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data1.out')
##data2 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data2.out')
data3 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data3.out')
data4 = np.loadtxt('RG_final_data/08_27Ts/cotangent_data/RG_cotangent_data4.out')
data5 = np.loadtxt('RG_final_data/09_250Ts/cotangent_data/RG_cotangent_data5.out')
##data6 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data6.out')
##data7 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data7.out')
##data8 = np.loadtxt('RG_final_data/06_20Ts_fixed/cotangent_data/RG_cotangent_data8.out')




numdata=1000
max_fillingfraction_handled = 0.55

x = plt.linspace(1e-20,max_fillingfraction_handled/(4*np.pi/3),numdata)                         # x-axis grid space (used for when plotting number density on x-axis)
xmin=(3/(4.0*np.pi))*(1e-20)*(1/(SW.R)**3)
xmax=(3/(4.0*np.pi))*(0.2)*(1/(SW.R)**3)
xff = plt.linspace(xmin,xmax,numdata)



  

##nL1 = [data1[i][1] for i in range(0,len(data1))]           # list of the left number density solutions obtained from fsolve
##nR1 = [data1[i][2] for i in range(0,len(data1))]           # list of the right number density solutions obtained from fsolve
##Tlist1 = [data1[i][0] for i in range(0,len(data1))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL1 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL1]           # converts left number density to filling fraction
##ffR1 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR1]           # converts right number density to filling fraction


##nL2 = [data2[i][1] for i in range(0,len(data2))]           # list of the left number density solutions obtained from fsolve
##nR2 = [data2[i][2] for i in range(0,len(data2))]           # list of the right number density solutions obtained from fsolve
##Tlist2 = [data2[i][0] for i in range(0,len(data2))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL2 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL2]           # converts left number density to filling fraction
##ffR2 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR2]           # converts right number density to filling fraction
##
##
nL1 = [data3[i][1] for i in range(0,len(data3))]           # list of the left number density solutions obtained from fsolve
nR1 = [data3[i][2] for i in range(0,len(data3))]           # list of the right number density solutions obtained from fsolve
Tlist1 = [data3[i][0] for i in range(0,len(data3))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL3 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL3]           # converts left number density to filling fraction
##ffR3 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR3]           # converts right number density to filling fraction


##nL1 = [data4[i][1] for i in range(0,len(data4))]           # list of the left number density solutions obtained from fsolve
##nR1 = [data4[i][2] for i in range(0,len(data4))]           # list of the right number density solutions obtained from fsolve
##Tlist1 = [data4[i][0] for i in range(0,len(data4))]        # list of the corresponding temperatures for which the above number densitites were found
##
##ffL4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL4]           # converts left number density to filling fraction
##ffR4 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR4]           # converts right number density to filling fraction
##
##
##nL1 = [data5[i][1] for i in range(0,len(data5))]           # list of the left number density solutions obtained from fsolve
##nR1 = [data5[i][2] for i in range(0,len(data5))]           # list of the right number density solutions obtained from fsolve
##Tlist1 = [data5[i][0] for i in range(0,len(data5))]        # list of the corresponding temperatures for which the above number densitites were found

##ffL5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nL5]           # converts left number density to filling fraction
##ffR5 = [i*((4*np.pi*(SW.R)**3)/3) for i in nR5]           # converts right number density to filling fraction










numdensity = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),numdata)

numdensity2 = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),numdata)

numdata2=250    #Number of numdensity data points
max_fillingfraction_handled2 = 0.15
ffs = plt.linspace(0.0001,max_fillingfraction_handled2,numdata2)
T = []

for Ts in Tlist1:
	T.append(Ts)




data = []

#for p in range(23,numdata2):
for p in range(45,len(T)):

        temp=T[p]
        print temp




        def fload():
                # Loads RG data and splits the file into two lists for number density and RG free energy
                # Also uses the interp1d function from scipy to interpolate between data points
                global finterp


                fname = 'RG_final_data/09_250Ts/data2/fit_T%.16g_f3.out'%temp
                #print fname," %.17g"%temp
                if(os.path.isfile(fname)==False):
                        fname = 'RG_final_data2/09_250Ts/data/fit_T%.16g_f3.out'%(temp+1e-16)
                if(os.path.isfile(fname)==False):
                        fname = 'RG_final_data2/09_250Ts/data/fit_T%.16g_f3.out'%(temp-1e-16)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%temp
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp+1e-15)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp-1e-15)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp+1e-14)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp-1e-14)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp+1e-13)
                if(os.path.isfile(fname)==False):
                        fname='RG_final_data2/09_250Ts/data/fit_T%.15g_f3.out'%(temp-1e-13)
                
                fdata=np.loadtxt(fname)
                #print(len(fdata))
                f = [fdata[i][1] for i in range(0,len(fdata))]
                #print(f[0])
                numdensity = [fdata[i][0] for i in range(0,len(fdata))]

                
                #print(numdensity[0])
                finterp = interp1d(numdensity,f,kind='cubic')



        def f_tot(n):
                #  Returns total free energy data with the missing a1SW(n) term that was not used in RG process
                return f_ext(n) + RG.a1SW(n)*n


        def f_ext(numdensity):
                # Fsolve might need to go outside of our minimum or maximum density that RG was used with, thus this function...
                # ... sets the values of the free energy outside those bounds.
                if numdensity > 0.0001 and numdensity < max_fillingfraction_handled/(4*np.pi/3):
                        return finterp(numdensity)
                return RG.fiterative(temp,numdensity,0)


        def dfext_dn(n):
                dn = 1e-6*n
                return (f_ext(n + dn) - f_ext(n - dn))/(2*dn)
                




        def frgp_ext(numdensity):
                # Fsolve might need to go outside of our minimum or maximum density that RG was used with, thus this function...
                # ... sets the values of the free energy outside those bounds. 
                if numdensity > 0.0001 and numdensity < max_fillingfraction_handled/(4*np.pi/3):
                        return frgpinterp(numdensity)
                if numdensity > max_fillingfraction_handled/(4*np.pi/3):
                        return frgpinterp(max_fillingfraction_handled/(4*np.pi/3))
                return frgpinterp(0.0001)
        

        def dfrgp_dn(n):
            ### Derivative of total free energy per volume w.r.t. n  (NOTE: at fixed volume and pressure, this is mu) ###
            dn = 1e-6*n
            return (frgp_ext(n + dn) - frgp_ext(n - dn))/(2*dn)










        def plotstuff():
                global rgp
                global frgpinterp
                numdensity3 = plt.linspace(0.0001,max_fillingfraction_handled/(4*np.pi/3),1000)
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
              

                yrgp = []
                dyrgp = []
                madyrgp = []
                for i in range(0,len(numdensity3)):
                        yrgp.append(frgp_ext(numdensity3[i]))
                        dyrgp.append(dfrgp_dn(numdensity3[i]))
                        if numdensity3[i] < 0.085 and numdensity3[i]>0.001:
                                madyrgp.append(dyrgp[i])
                madyrgp = sum(madyrgp)/len(madyrgp)
                mu=SW.numH_dftot_dn(T[p],numdensity3)
                mu2 = SW.numH_dftot_dn(T[p],nR1[p])
                mu3 = dfrgp_dn(nL1[p])
                mu4 = dfrgp_dn(nR1[p])
                #print madyrgp
                cotangentGuess = []
                for i in range(0,len(numdensity3)):
                        cotangentGuess.append(yrgp[i]-madyrgp*numdensity3[i])

                gfe=[]
                for i in range(0,len(numdensity3)):
                        gfe.append(yrgp[i]-mu4*numdensity3[i])


##                gfergp=[]
##                for i in range(0,len(numdensity3)):
##                        gfergp.append(yrgp[i]-mu3*numdensity3[i])



                swgfe=[]
                for i in range(0,len(numdensity3)):
                        swgfe.append(ysw[i]-mu2*numdensity3[i])
                        

##                plt.plot(numdensity3,cotangentGuess)
                plt.plot(numdensity3,swgfe,color='b',linewidth=2)
                plt.plot(numdensity3,gfe,color='#f36118',linewidth=2)
##                plt.plot(numdensity3,gfergp,color='c',linewidth=2)
##                plt.plot(numdensity3,gfergp,color='g',linewidth=2)
                plt.plot(nL1[p],gfe[p],'ko')
                plt.plot(nR1[p],gfe[p],'ko')
##                plt.axhline(SW.phi(T[p],nR1[p],nR1[p]),color='c',linewidth=2)
                plt.axhline(gfe[p],color='r',linewidth=2)

                plt.xlim(0,nR1[p]+0.01)
                plt.ylim(-.010,0.0075)
##                plt.show()
                plt.savefig('meeting/18july2016/f03/gfe_T%.16g.png'%T[p])
                plt.close()
              



        plotstuff()





        


        









