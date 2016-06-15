from __future__ import division
import os, sys, re
import numpy as np
from numpy import exp,log
import matplotlib
import pylab as plt
from itertools import cycle
useCarnahanStarling=True

figDir=os.getcwd()

sys.path.append(figDir+"/../../thesis-roth/figs/")
import SW

sys.path.append(figDir+"/../data")
import checkLiquidVapor as vapor
import thermoStats
import listAll
listAll.baseDir+="/../data"
listAll.findFolders()

def SW_Fexc(T,n):
	#return SW.fhs(T,n)
	#return SW.fdisp(T,n)
	return SW.fdisp(T,n)+SW.fhs(T,n)

def SW_Finfinity(n):
	T=1.0e8
	return SW.fhs(T,n)/T

def avgFcost(stats,Tmin,Tmax):
	T=Tmin
	deltaT=(Tmax-Tmin)/1000
	count=0
	totalCost=0
	while T<Tmax:
		totalCost+=Fcost(stats,T)
		T+=deltaT
		count+=1
	return totalCost/count
def FcostFixedFF(stats):
	ffs=[0.05,0.15,0.25,0.35,0.45]
	nPairs=[]
	for ff in ffs:#first find the nPairs corresponding to each FF############
		nLeft=int(ff*(stats[0].L**3.0)/(4/3.0*np.pi))#assumes R=1.0
		nRight=int(ff*(stats[0].L**3.0)/(4/3.0*np.pi)+1.0)
		nPair=[]
		for i in range(0,len(stats)):
			if stats[i].N==nLeft:
				nPair.append(i)
				#
				break
		for i in range(0,len(stats)):
			if stats[i].N==nRight:
				nPair.append(i)
				#
				break
		if len(nPair)<2: return -1.0
		nPairs.append(nPair)
	if len(nPairs)==0: return -1.0
	diffSquared=0.0
	weightTotal=0
	Ts=plt.linspace(0.6,1.3,1000)
	for i in range(0,len(nPairs)):#now find the diffSquared by interpolation of the pairs############
		nMid=ffs[i]*(stats[0].L**3.0)/(4/3.0*np.pi) #what we want
		nLeft=int(ffs[i]*(stats[0].L**3.0)/(4/3.0*np.pi))#assumes R=1.0
		density=nMid/stats[0].L**3.0
		for T in Ts:
			Fleft=stats[nPairs[i][0]].findFexc(T)
			Fright=stats[nPairs[i][1]].findFexc(T)
			Fmid=Fleft+(Fright-Fleft)/(1.0)*(nMid-nLeft) #F0+(rise/run)*dn
			diffSquared+=(Fmid-SW_Fexc(T,density))**2.0/density**2.0
			weightTotal+=1
	return (diffSquared/weightTotal)**0.5

		
def Fcost(stats,T):
	#free energy cost function
	#now uses RMS
	diffSquared=0.0
	weightTotal=0.0
	for thermoStat in stats:
		density=thermoStat.N/(thermoStat.L**3.0)
		diffSquared+=((thermoStat.findFexc(T)-SW_Fexc(T,density))**2.0)/density**2.0
		weightTotal+=1
	return (diffSquared/weightTotal)

dataColors=['ro','go','bo','rx','gx','bx']
saftColors=['r','g','b','r--','b--','g--']

#iterates over base folders and ith orders
#finally loads the N atoms in a box data
#for each set, generates some plots near optimal temp (found in 1.0->1.3
LforAvgCost=[]
avgCost=[]
for m in range(0,len(listAll.validFoldersss)):#iterates over base folders
	for ith in range(0,len(listAll.validFoldersss[m])):#iterates over ith orders
		stats=[]	#holds thermoStats
		for nj in range(0,len(listAll.validFoldersss[m][ith])):#loads N_j data
			fbase=listAll.validFoldersss[m][ith][nj]
			if False==os.path.isfile(fbase+"lv-data-dos.dat"): continue
			if useCarnahanStarling==True:
				print fbase
				L=float(fbase.split("-L")[1].split("/")[0])
				density=int(fbase.split("N")[1].split("/")[0])/L**3.0
				temp=thermoStats.ThermoStats(fbase=fbase,SatInfinity=-SW_Finfinity(density)*(L**3.0))
				
			else:
				if False==os.path.isfile(fbase+"absolute/Sexc.dat"): continue
				temp=thermoStats.ThermoStats(fbase=fbase)
			#below taken out since final data does not inlcude the .out file
			#left in to show I do like to check if data is ready...
			#suggestion: don't generate Sexc.dat until liquidVapor is done
			#if False==vapor.check(fbase=fbase): continue
			
			if len(temp.dos)<1: continue
			stats.append(temp)
			print temp.L
		if len(stats)==0: continue
		#data now loaded, PLOT!
		
		####------ plots free energy around the optimal temp ---########
		'''maxFF=0.5
		density=plt.linspace(0.01,maxFF*3/4.0/np.pi,1000)
		ff=density*4/3.0*np.pi
		plt.figure()
		ffData=[]
		FexcData=[]
		UexcData=[]
		SexcData=[]
		for j in range(0,len(stats)):
			ffData.append(stats[j].N*4*np.pi/3.0/((stats[j].L**3.0)))
		Ts=plt.linspace(1.0,1.3,1000)
		Toptimal=(Ts[0]+Ts[-1])/2.0
		SW.lambdaSW=stats[0].ww
		print SW.lambdaSW
		FcostOptimal=Fcost(stats,Toptimal)
		for T in Ts:
			fcost=Fcost(stats,T)
			if fcost<FcostOptimal:
				Toptimal=T
				FcostOptimal=fcost
		temperatures=[]
		if 0.8*Toptimal>0.55:
			temperatures.append(0.8*Toptimal)
		temperatures.append(Toptimal)
		temperatures.append(1.2*Toptimal)
		for j in range(0,len(temperatures)):
			SWdata=[]
			for jj in range(0,len(density)):
				SWdata.append(SW_Fexc(temperatures[j],density[jj])/density[jj]) #convert to per N
			plt.plot(ff,SWdata,saftColors[j%len(saftColors)])
			yData=[]
			for k in range(0,len(stats)):
				yData.append(stats[k].findFexc(temperatures[j])*(stats[k].L**3.0)/stats[k].N)
			plt.plot(ffData,yData,dataColors[j%len(dataColors)],label="T=%f"%temperatures[j])
		title="FexcVsff-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("filling fraction")
		plt.ylabel("$F_{exc}$")
		plt.title(title)
		plt.legend(loc='best')
		plt.savefig(title+".pdf")
		plt.close()'''
		
		##################-------------COST GRAPH----------#############
		#we stats[n_i], but we only want FF=[0.05,0.15,0.25,0.35,0.45]
		'''Ts=plt.linspace(1.0,1.3,1000)
		Fcosts=[]
		for T in Ts:
			Fcosts.append(Fcost(stats,T))
		plt.figure()
		plt.plot(Ts,Fcosts,'r',label='raw Data')
		edgeX0,edgeX1,edgeY0,edgeY1=plt.axis()
		#plt.axis((edgeX0-0.2,edgeX1-0.2,edgeY0-0.2,edgeY1+0.1))
		plt.legend(loc='best')
		title="FcostVsT-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("Temperature")
		plt.ylabel("RMS(Fexc-SaftF)")#"$\\frac{\\Sigma((Fexc-SaftF)/\\rho)^2}{\\Sigma\\rho^2}$")
		plt.title(title)
		plt.legend(loc='best')
		plt.savefig(title+".pdf")
		plt.close()'''
		
		if abs(stats[0].ww-1.5)<1e-4:
			currentCost=FcostFixedFF(stats)
			if currentCost>=0:
				LforAvgCost.append(stats[0].L)
				avgCost.append(currentCost)

###############--------------all ww=1.5  T=infinity  dataF-saftF------##
colors=['ro','go','bo','rx','gx','bx']
plt.figure()
plt.plot(LforAvgCost,avgCost,'ro')
plt.xlabel("length of box")
plt.ylabel("average RMS(Fexc-SaftF) for T in (0.7,1.5)")
plt.savefig("costVsBoxSize.png")
plt.close()
