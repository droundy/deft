from __future__ import division
import os, sys, re
import numpy as np
from numpy import exp,log
import matplotlib
import pylab as plt
from itertools import cycle

figDir=os.getcwd()

sys.path.append(figDir+"/../../thesis-roth/figs/")
import SW

sys.path.append(figDir+"/../data")
import checkLiquidVapor as vapor
import thermoStats
import listAll
listAll.baseDir+="/../data"
#print listAll.baseDir
listAll.findFolders()

def SW_Fexc(T,n):
	#return SW.fhs(T,n)
	#return SW.fdisp(T,n)
	return SW.fdisp(T,n)+SW.fhs(T,n)

def SW_Finfinity(n):
	T=1.0e8
	return SW.fhs(T,n)/T

def Fcost(stats,T):
	#free energy cost function
	diffSquared=0.0
	weightTotal=0.0
	for thermoStat in stats:
		density=thermoStat.N/(thermoStat.L**3.0)
		diffSquared+=((thermoStat.findFexc(T)-SW_Fexc(T,density))**2.0)*(density**2.0)
		weightTotal+=density**2.0
	return diffSquared/weightTotal

def FcostGuessFinfinity(stats,T):
	#free energy cost function
	#uses saft to find Sexc at T=infinity
	diffSquared=0.0
	weightTotal=0.0
	for thermoStat in stats:
		density=thermoStat.N/(thermoStat.L**3.0)
		sTemp=thermoStat.SatInfinity
		thermoStat.SatInfinity=-(SW_Finfinity(density))*(thermoStat.L**3.0)   #-0.005/(0.35*3.0/4.0/np.pi)*density
		diffSquared+=((thermoStat.findFexc(T)-SW_Fexc(T,density))**2.0)*(density**2.0)
		thermoStat.SatInfinity=sTemp
		weightTotal+=density**2.0
	return diffSquared/weightTotal

dataColors=['ro','go','bo','rx','gx','bx']
saftColors=['r','g','b','r--','b--','g--']

#iterates over base folders and ith orders
#finally loads the N atoms in a box data
#for each set, generates some plots near optimal temp (found in 1.0->1.3
allFinfinity=[] #used to print all the -Sexc
allFinfinityPercent=[] #difference between theory and data
allCellLength=[] #cell length for above data (otherwise they get thrown away)
for m in range(0,len(listAll.validFoldersss)):#iterates over base folders
	for ith in range(0,len(listAll.validFoldersss[m])):#iterates over ith orders
		stats=[]	#holds thermoStats
		for nj in range(0,len(listAll.validFoldersss[m][ith])):#loads N_j data
			fbase=listAll.validFoldersss[m][ith][nj]
			if False==os.path.isfile(fbase+"absolute/Sexc.dat"): continue
			#below taken out since final data does not inlcude the .out file
			#left in to show I do like to check if data is ready...
			#suggestion: don't generate Sexc.dat until liquidVapor is done
			#if False==vapor.check(fbase=fbase): continue
			if False==os.path.isfile(fbase+"lv-data-dos.dat"): continue
			temp=thermoStats.ThermoStats(fbase=fbase)
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
			plt.plot(ff,SW_Fexc(temperatures[j],density),saftColors[j%len(saftColors)])
			yData=[]
			for k in range(0,len(stats)):
				yData.append(stats[k].findFexc(temperatures[j]))
			plt.plot(ffData,yData,dataColors[j%len(dataColors)],label="T=%f"%temperatures[j])
		title="FexcVsff-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("filling fraction")
		plt.ylabel("$F_{exc}$")
		plt.title(title)
		plt.legend(loc='best')
		plt.savefig(title+".pdf")
		plt.close()'''
		
		# plots both free energy and mod free energy using Sexc(T=infinity) from saft ###
		maxFF=0.5
		density=plt.linspace(0.01,maxFF*3/4.0/np.pi,1000)
		ff=density*4/3.0*np.pi
		plt.figure(figsize=(16,6),dpi=300)
		plt.subplot(121)
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
			plt.plot(ff,SW_Fexc(temperatures[j],density),saftColors[j%len(saftColors)])
			yData=[]
			for k in range(0,len(stats)):
				yData.append(stats[k].findFexc(temperatures[j]))
			plt.plot(ffData,yData,dataColors[j%len(dataColors)],label="T=%f"%temperatures[j])
		title="FexcVsff-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("filling fraction")
		plt.ylabel("$F_{exc}$")
		plt.title(title)
		plt.legend(loc='best')
		
		plt.subplot(122)
		#Toptimal=(Ts[0]+Ts[-1])/2.0
		SW.lambdaSW=stats[0].ww
		print SW.lambdaSW
		'''FcostOptimal=FcostGuessFinfinity(stats,Toptimal)
		for T in Ts:
			fcost=FcostGuessFinfinity(stats,T)
			if fcost<FcostOptimal:
				Toptimal=T
				FcostOptimal=fcost
		temperatures=[]
		if 0.8*Toptimal>0.55:
			temperatures.append(0.8*Toptimal)
		temperatures.append(Toptimal)
		temperatures.append(1.2*Toptimal)'''
		for j in range(0,len(temperatures)):
			plt.plot(ff,SW_Fexc(temperatures[j],density),saftColors[j%len(saftColors)])
			yData=[]
			for k in range(0,len(stats)):
				sTemp=stats[k].SatInfinity
				stats[k].SatInfinity=-(SW_Finfinity(stats[k].N/stats[k].L**3.0))*(stats[k].L**3.0)
				yData.append(stats[k].findFexc(temperatures[j]))
				stats[k].SatInfinity=sTemp
			plt.plot(ffData,yData,dataColors[j%len(dataColors)],label="T=%f"%temperatures[j])
		title="modified using SW to find S at infinity"#"FexcVsff-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("filling fraction")
		plt.ylabel("$F_{exc}$")
		plt.title(title)
		plt.legend(loc='best')
		
		title="FexcVsff-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.savefig(title+".pdf",dpi=300)
		plt.close()
		
		##################-------------COST GRAPH----------#############
		Ts=plt.linspace(1.0,1.3,1000)
		Fcosts=[]
		for T in Ts:
			Fcosts.append(Fcost(stats,T))
		plt.figure()
		plt.plot(Ts,Fcosts,'r',label='raw Data')
		Fcosts=[]
		for T in Ts:
			Fcosts.append(FcostGuessFinfinity(stats,T))
		plt.plot(Ts,Fcosts,'b',label='guess Sinfinity')
		edgeX0,edgeX1,edgeY0,edgeY1=plt.axis()
		#plt.axis((edgeX0-0.2,edgeX1-0.2,edgeY0-0.2,edgeY1+0.1))
		plt.legend(loc='best')
		title="FcostVsT-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("Temperature")
		plt.ylabel("$\\frac{\\Sigma((Fexc-SaftF)*\\rho)^2}{\\Sigma\\rho^2}$")
		plt.title(title)
		plt.legend(loc='best')
		plt.savefig(title+".pdf")
		plt.close()
		
		#################------Free energy at infinity graphs--#########
		plt.figure()
		freeEnergySaft=[]
		maxFF=0.5
		density=plt.linspace(0.01,maxFF*3/4.0/np.pi,1000)
		ff=density*4/3.0*np.pi
		for i in range(0,len(density)):
			freeEnergySaft.append(SW_Finfinity(density[i]))
		plt.plot(ff,freeEnergySaft,'r')
		dataFreeEnergy=[]
		ffData=[]
		for i in range(0,len(stats)):
			dataFreeEnergy.append(-stats[i].SatInfinity/(stats[i].L**3.0))
			ffData.append(stats[i].N*4*np.pi/3.0/((stats[i].L**3.0)))
		plt.plot(ffData,dataFreeEnergy,'go')
		title="Finfinity-"+fbase.split("data/")[1].split("/N")[0].replace('/','-')
		plt.xlabel("filling fraction")
		plt.ylabel("$F_{exc}/T$  for  T->infinity")
		plt.savefig(title+".pdf")
		plt.close()
		#prepare combination plot comparison that includes ALL ww=1.5
		if abs(stats[0].ww-1.5)<1e-4:
			allCellLength.append(stats[0].L)
			for i in range(0,len(dataFreeEnergy)):
				dataFreeEnergy[i]-=SW_Finfinity(stats[i].N/(stats[i].L**3.0))
			#allFinfinity.append([ffData,dataFreeEnergy])
			allFinfinity.append([[],[]])
			for i in range(0,len(dataFreeEnergy)):
				allFinfinity[-1][0].append(ffData[i])
				allFinfinity[-1][1].append(dataFreeEnergy[i])
			
			for i in range(0,len(dataFreeEnergy)):
				dataFreeEnergy[i]=abs(dataFreeEnergy[i]/stats[i].SatInfinity)*100.0
			#allFinfinityPercent.append([ffData,dataFreeEnergy])
			allFinfinityPercent.append([[],[]])
			for i in range(0,len(dataFreeEnergy)):
				allFinfinityPercent[-1][0].append(ffData[i])
				allFinfinityPercent[-1][1].append(dataFreeEnergy[i])
			
###############--------------all ww=1.5  T=infinity  dataF-saftF------##
colors=['ro','go','bo','rx','gx','bx']
plt.figure()
size=15.0
for i in range(0,len(allFinfinity)):
	plt.plot(allFinfinity[i][0],allFinfinity[i][1],colors[i%len(colors)],label="L=%f"%allCellLength[i])
plt.xlabel("filling fraction")
plt.ylabel("dataF-saftF")
plt.legend(loc='best')
plt.savefig("deltaF.pdf")
plt.close()

colors=['ro','go','bo','rx','gx','bx']
plt.figure()
size=15.0
for i in range(0,len(allFinfinity)):#ms=size*(0.7)**i
	plt.plot(allFinfinityPercent[i][0],allFinfinityPercent[i][1],colors[i%len(colors)],label="L=%f"%allCellLength[i])
plt.xlabel("filling fraction")
plt.ylabel("(dataF-saftF)/dataF*100.0")
plt.legend(loc='best')
plt.savefig("deltaFpercent.pdf")
plt.close()

print SW.fhs(0,0.1)
print SW.fhs(0,0.2)
