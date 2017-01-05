#!/usr/bin/env python2
from __future__ import division
import os, sys, re
import numpy as np
from numpy import exp,log
import matplotlib
import pylab as plt
from itertools import cycle

class ThermoStats:
    def __init__(self,fbase=None,SatInfinity=None):
        if fbase==None: return
        if fbase[-1]!='/': fbase+='/'
        self.ww=None
        self.L=None
        self.N=None
        self.SatInfinity=None
        if SatInfinity!=None:
			self.SatInfinity=SatInfinity
        else:
			if os.path.isfile(fbase+"absolute/Sexc.dat")==False: return
			self.SatInfinity=float(np.loadtxt(fbase+"absolute/Sexc.dat"))
        dataList=np.loadtxt(fbase+'lv-data-dos.dat')
        dos=dataList[:,1]
        energy=dataList[:,0]
        #print fbase
        if len(dos)>1 and dos[0]==dos[1]:
	        while len(dos)>1 and dos[0]==dos[1]:
	            dos=np.delete(dos,0,axis=0)
	            energy=np.delete(energy,0,axis=0)
	        dos=np.delete(dos,0,axis=0)
	        energy=np.delete(energy,0,axis=0)
        if len(dos)>1 and dos[-1]==dos[-2]:
	        while len(dos)>1 and dos[-1]==dos[-2]:
	            dos=np.delete(dos,-1,axis=0)
	            energy=np.delete(energy,-1,axis=0)
	        dos=np.delete(dos,-1,axis=0)
	        energy=np.delete(energy,-1,axis=0)
        self.dos=dos
        self.energy=energy
        ww=None
        L=None
        N=None
        with open(fbase+'lv-data-dos.dat','r') as file:
            for line in file:
                if 'N:' in line: N=int(line.split()[-1])
                if 'cell dimensions' in line: L=float(line.split()[-1].split(")")[0])
                if 'well_width' in line: ww=float(line.split()[-1])
        if L==None: L=-1.0; print "did not find L"
        if ww==None: ww=-1.0; print "did not find ww"
        if N==None: N=-1.0; print "did not find N"
        self.ww=ww
        self.ff=4/3.0*np.pi*N/(L**3.0)
        self.L=L
        self.N=N
        return
    def findSexc(self,T):
        """ returns Uexc(n,T)/volume/T - (Fexc(n,T)/volume/T - Fexc(n,T=infinity)/volume/(T=infinity))
                    = Uexc(n,T)/volume/T - Fexc(n,T)/volume/T - Sexc(n,T=infinity)/volume
                    = Sexc(n,T)/volume - Sexc(n,T=infinity)/volume
                    = Sdisp(n,T)/volume - Sdisp(n,T=infinity)/volume """
        return (self.findUexc(T)-self.findFexc(T))/T
    
    def findUexc(self,T):
        """ return the excess internal energy (for each temperature
            in T), per unit volume, in
            units where the radius of the sphere is 1. """
        energy=self.energy
        dos=self.dos
        z=dos+energy/T
        energyWeights=np.log(energy)+dos+energy/T
        zMax=np.max(z)
        weightMax=np.max(energyWeights)
        energyWeights=np.sum(exp(energyWeights-weightMax))
        z=np.sum(np.exp(z-zMax))
        return -np.exp(weightMax-zMax)*energyWeights/z/(self.L**3.0)
    
    def findFexc(self,T):
        """ returns Fexc(n,T)/volume - T*Fexc(n,T=infinity)/volume/(T=infinity)
                    = Fexc(n,T)/volume + T*Sexc(n,T=infinity)/volume """
        energy=self.energy
        dos=self.dos
        z=dos+energy/T
        zMax=z.max()
        dosMax=np.max(dos)
        zInfinity=np.log(np.sum(np.exp(dos-dosMax)))+dosMax
        return (-T*(log(np.sum(exp(z-zMax)))+zMax-zInfinity)-T*self.SatInfinity)/(self.L**3.0)
        
    def generateData(self,T):
        Uexc=np.zeros(len(T))
        Sexc=np.zeros(len(T))
        Fexc=np.zeros(len(T))
        for i in range(0,len(T)):
            Uexc[i]=self.findUexc(T[i])
            Sexc[i]=self.findSexc(T[i])
            Fexc[i]=self.findFexc(T[i])
        self.Uexc=Uexc
        self.Sexc=Sexc
        self.Fexc=Fexc
        return


if __name__=='__main__':
    print "hi"
    

