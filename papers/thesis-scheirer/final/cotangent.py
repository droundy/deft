import scipy as sp
from scipy.optimize import fsolve
import pylab as plt
import matplotlib
import SW
import numpy as np


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
temp = plt.linspace(0.3,0.3,1)

### non-uniform temperature linspace (used to get higher resolution near critical temperature where fsolve starts to break down) ###
#temp = np.concatenate((plt.linspace(0.3,1.0,300),plt.linspace(1.0+0.7/300,1.3,600),plt.linspace(1.3+0.3/600,1.328,1000),plt.linspace(1.328+0.028/1000,1.33,30000)),axis=0)

### Initial guesses for fsovle ###
### (NOTE:  need to "play" with these initally to get fsolve to work... ###
### ...This may requires plotting your free energy per volume vs number density first to come up with some initial guesses. ###
#Lguess = 1e-9  # Initial left "n" guess (good for 0.6)
Lguess = 1e-12  # Initial left "n" guess (good for 0.3)
Rguess = 0.2    # Initial right "n" guess
#                                                                                             #
#                                                                                             #
############################### END INITIALIATION #############################################





############################### START fsolve WORK #############################################
#                                                                                             #
#                                                                                             #
 ###
def g(x):
        ### just a general function which includes conditions for finding the cotangent
	x1=x[0]                                 # left guess "n"
	x2=x[1]                                 # right guess "n"
	y1=SW.ftot(T,x[0])                      # f(x1)
	y2=SW.ftot(T,x[1])                      # f(x2)
	dy1=SW.numH_dftot_dn(T,x[0])            # df_dx(x1)
	dy2=SW.numH_dftot_dn(T,x[1])            # df_dx(x2)

                             
	out=[(dy1-dy2)]                         # Condition 1: df_dx(x1) = df_dx(x2)
	out.append(dy1-((y2-y1)/(x2-x1)))       # Condition 2: df_dx(x1) = slope between the two positions
	return out

data=[]
count=0
for T in temp:
	count+=1
	if count%10==0:
		print '%f    \r'%T,
	sol = fsolve(g,[Lguess,Rguess])         # Magic happens here
	data.append([T,sol[0],sol[1]])          # updating a list to keep track of the temperature T and corresponding left and right number densities found by fsolve
	Lguess=sol[0]                           # Resets the next left guess number density for fsolve as the previous left solution
	Rguess=sol[1]                           # Resets the next right guess number density for fsolve as the previous right solution
print(sol)
print(temp[-1])
#                                                                                             #
#                                                                                             #
########################### END fsolve WORK ###################################################	





########################### START PLOTTING STUFF ##############################################        
#                                                                                             #
#                                                                                             #
x = plt.linspace(1e-20,.2,4000)                         # x-axis grid space (used for when plotting number density on x-axis)
xmin=(3/(4.0*np.pi))*(1e-20)*(1/(SW.R)**3)
xmax=(3/(4.0*np.pi))*(0.2)*(1/(SW.R)**3)
xff = plt.linspace(xmin,xmax,4000)


nL = [data[i][1] for i in range(0,len(data))]           # list of the left number density solutions obtained from fsolve
nR = [data[i][2] for i in range(0,len(data))]           # list of the right number density solutions obtained from fsolve
Tlist = [data[i][0] for i in range(0,len(data))]        # list of the corresponding temperatures for which the above number densitites were found

ffL = [i*((4*np.pi*(SW.R)**3)/3) for i in nL]           # converts left number density to filling fraction
ffR = [i*((4*np.pi*(SW.R)**3)/3) for i in nR]           # converts right number density to filling fraction



def cotangent_t0():
        ### Total free energy per volume  VS  n  (with cotangent line shown) for the first temperature  ###
        plt.figure()
        plt.title('Total free energy per volume  VS  n  @ T=%0.4f'%Tlist[-1])
        plt.ylabel('Total free energy per volume')
        plt.xlabel('Number density (n)')
        plt.plot(x,SW.ftot(Tlist[-1],x),color='#f36118',linewidth=1)
        plt.plot(x, SW.numH_dftot_dn(Tlist[-1],nR[-1])*(x-nR[-1])+SW.ftot(Tlist[-1],nR[-1]),color='#00c0c0',linewidth=1)
        plt.plot(nL[-1],SW.ftot(Tlist[-1],nL[-1]),'ko')
        plt.plot(nR[-1],SW.ftot(Tlist[0],nR[-1]),'ko')
        #plt.legend(loc='best')
        plt.savefig('cotangent.pdf')
        plt.show()
        #plt.close("all")



def cotangent_tloop():
        ### Total free energy per volume  VS  n  (with cotangent line shown) for multiple temperatures, can use to make a gif later ###
        count = 0 
        for i in range(0,len(Tlist)):
                if (i%20==0):
                        plt.figure()
                        plt.title('Total free energy per volume  VS  n  @ T=%0.4f'%Tlist[i])
                        plt.ylabel('Total free energy per volume')
                        plt.xlabel('Number density (n)')
                        plt.ylim(-0.8,1)
                        plt.xlim(0,.18)
                        plt.plot(x,SW.ftot(Tlist[i],x))
                        plt.plot(x, SW.numH_dftot_dn(Tlist[i],nR[i])*(x-nR[i])+SW.ftot(Tlist[i],nR[i]))
                        plt.plot(nL[i],SW.ftot(Tlist[i],nL[i]),'ko')
                        plt.plot(nR[i],SW.ftot(Tlist[i],nR[i]),'ro')
                        plt.savefig('cotangent_loop/cotangent%03d'%count)
                        count += 1



def liq_vap_co_tvsff():
        ### Vapor - Liquid coexistence curve for temperature  VS  filling fraction ###
        plt.figure()
        plt.title('Vapor - Liquid coexistence')
        plt.ylabel('T')
        plt.xlabel('ff')
        plt.plot(ffL,Tlist,'m')
        plt.plot(ffR,Tlist,'c')
        plt.savefig('liqVapCo_Tvsff.pdf')
        #plt.axhline(1.329,color='k')                           
        #plt.axhline(1.33,color='r')
        #plt.show()


def liq_vap_co_tvsn():
        ### Vapor - Liquid coexistence curve for temperature  VS  filling fraction ###
        plt.figure()
        plt.title('Vapor - Liquid coexistence')
        plt.ylabel('T')
        plt.xlabel('ff')
        plt.plot(nL,Tlist,'m')
        plt.plot(nR,Tlist,'c')
        plt.savefig('liqVapCo_Tvsn.pdf')
        #plt.axhline(1.329,color='k')                           
        #plt.axhline(1.33,color='r')
        #plt.show()



def cotangent_fsolve_breakdown():
        ### Plots free energy per volume scaled by the tangent found from fsolve (used to determine approximately when fsolve stops working) ###
        previousSize = nR[0]-nL[0]
        count = 0

        for i in range(0,len(nL)):
                if (nR[i]-nL[i]) > 0.001: continue
        	#if (nR[i]-nL[i]) < 0.0000001: break
                if True:
                #if previousSize - (nR[i]-nL[i]) > 0.0000001:
                        x = plt.linspace(nL[i],nR[i],4000)
                        previousSize = nR[i]-nL[i]
                        datatangent=SW.numH_dftot_dn(Tlist[i],nR[i])*(x-nR[i])+SW.ftot(Tlist[i],nR[i])
                        #Tprev = Tlist[i]
                        plt.figure()
                        plt.title('ftot scaled  VS  n  @ T=%0.6f'%Tlist[i])
                        plt.ylabel('Total free energy per volume')
                        plt.xlabel('n')		
                        plt.plot(x,SW.ftot(Tlist[i],x)-datatangent,'b')
                        plt.plot(x,x-x,'g')
                        plt.plot(nL[i],0,'ko')
                        plt.plot(nR[i],0,'ro')
                        #plt.legend(loc='best')
                        plt.savefig('cotangent-graphs/relativeH/graph%04d'%count)
                        count += 1
                        #print(nR[i]-nL[i])
                        plt.close("all")

                        

def liq_vap_co_pvsT():
        ### Vapor - Liquid coexistence curve for pressure vs temperature ###
        pL=[]
        pR=[]
        pdiff=[]
        for i in range(0,len(nL)):
                pL.append(SW.findP(Tlist[i],nL[i]))
                pR.append(SW.findP(Tlist[i],nR[i]))
                pdiff.append(pL[i]-pR[i])
        plt.figure()
        plt.title('Vapor - Liquid coexistence')
        plt.ylabel('Pressure')
        plt.xlabel('Temperature')
        plt.plot(Tlist,pL,color='#f36118',linewidth=5)
        plt.plot(Tlist,pR,'c',linewidth=1)
        plt.savefig('liqVapCo_pvsT.pdf')
        
        plt.figure()
        plt.title('Pressure check')
        plt.ylabel('P')
        plt.xlabel('T')
        plt.plot(Tlist,pdiff,'r')
        plt.savefig('liqVapCo_pvsT_diff.pdf')	



#################################################
#       CALL THE PLOTS YOU WANT TO PLOT         #
#                                               #
#cotangent_t0()
#cotangent_tloop()
#liq_vap_co_tvsff()
#liq_vap_co_tvsn()
#cotangent_fsolve_breakdown()
#liq_vap_co_pvsT()
#                                               #
#                                               #
#################################################

#np.savetxt('figs/snaft2.out',data)

#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
