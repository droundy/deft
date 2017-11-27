#!/usr/bin/env python2
from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import argparse
########################################################################
## This script reads data from files and plots for new-soft
########################################################################
#Global Constants
x = 0
y = 1
z = 2

# Default data sets to pull from. 
rhoDefault = [0.7,0.85,0.89,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,
        1.0,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,
        1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24]
TDefault = [2.0]

def plotPressure(rho,T):
	plt.figure()
	for density in rho:
		pressArray = []
		for temp in T:
			pressure = np.loadtxt('data2/ff-'+str(density)+\
                '_temp-'+str(temp)+'-press.dat')
                
                    
                pressArray.append(float(pressure))
                plt.plot(density,pressArray,'k.')
        plt.title('Pressure-Density at T*: '+str(temp))
        plt.ylabel('Pressure p*')
        plt.xlabel(r'Density $\rho*$')
        plt.grid()

def plotPositions(rho,T):
	for density in rho:
		for temp in T:
			spheres = np.loadtxt('data2/ff-'+str(density)+\
                '_temp-'+str(temp)+'-pos.dat')
			fig = plt.figure()
			ax = fig.add_subplot(111, projection = '3d')
			ax.scatter(spheres[:,0],spheres[:,1],spheres[:,2],\
                c = 'r', marker ='o')
			ax.set_title('Positions at Temp: '+ str(temp)+\
                ' and Density: ' + str(density))
			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')


def plotRadialDF(rho,T):
	for density in rho:
		for temp in T:
			radialData  = np.loadtxt('data2/ff-'+str(density)+\
                '_temp-'+str(temp)+'-radial.dat')
			plt.figure()
			plt.plot(radialData[:,0],radialData[:,1])#/(radialData[:,0]**2))
			plt.title('Radial Distribution Function.\n'\
				'Temp: '+str(temp)+' and Density: ' +str(density))
			plt.xlabel('Radial Distance (r)')
			plt.ylabel('Nearby Spheres')
                plt.grid()

def energy(rho,T):
    for temp in T:
        plt.figure()
        pressArray = []
        volumeArray = []
        for density in rho:
            fname = 'data2/ff-%s_temp-%s-press.dat'%(density,temp)
            pressArray.append(float(np.loadtxt(fname)))
            Nballs = -1
            with open(fname) as f:
                for l in f.readlines():
                    if l[:10] == '# volume: ':
                        volumeArray.append(float(l[10:]))
                        break
            with open(fname) as f:
                for l in f.readlines():
                    if l[:5] == '# N: ':
                        Nballs = int(l[5:])
                        break
        V = np.flip(volumeArray,0)
        p = np.flip(pressArray,0)
        plt.plot(V,p,'k.')
        plt.xlabel('Volume V')
        plt.ylabel('Pressure p*')
        plt.title('Pressure v. Volume At Constant Temp: %s'%temp)
        plt.grid()

        F = np.zeros(len(V)-1)
        pn,v = np.zeros(len(F)),np.zeros(len(F))
        F[0] = -(p[1]+p[0])*(V[1]-V[0])/(2*Nballs)
        v[0] = (V[1]+V[0])/(2*Nballs)
        pn[0] = (p[1]+p[0])/2
        for n in range(1,len(F)):
            F[n] += F[n-1]-(p[n+1]+p[n])*(V[n+1]-V[n])/(2*Nballs)
            v[n] = (V[n+1]+V[n])/(2*Nballs)
            pn[n] = (p[n+1]+p[n])/2

        plt.figure()
        plt.plot(v,F,'-.')
        plt.title(r'$f - f_0$')
        plt.xlabel('reduced volume v')
        plt.ylabel('reduced helmholtz')
        plt.grid()
        
        plt.figure()
        plt.plot(pn,F + pn*v - 2.8/3*pn,'k-')
        plt.xlabel('pressure p*')
        plt.ylabel('g(v)')
        plt.title(r' Attempted Reduced Gibbs with weird shift $g = (f-f_0) + pv$')
        plt.grid()
        
        rho = np.array(rho)
        mu = -22.0
        plt.figure()
        plt.plot(rho[:-1],F - mu*rho[:-1],'k-.')
        plt.xlabel('rho')
        plt.ylabel('g(v)')
        plt.title(r' Attempted grand with weird shift $g = (f-f_0) + pv$')
        plt.grid()
        
        
        
        
def plotStructureFactor(rho,T):
    
    for temp in T:
        for density in rho:
            fig, ax = plt.subplots()
            fname = 'data2/ff-%s_temp-%s-struc.dat' % (density, temp)
            structureData = np.loadtxt(fname)
            Nballs = -1;
            with open(fname) as f:
                for l in f.readlines():
                    if l[:5] == '# N: ':
                        Nballs = int(l[5:])
                        break
            dk = 1/40
            kvals = dk*np.arange(0.0, len(structureData), 1)
            kx, ky = np.meshgrid(kvals, kvals)
            resolution = 2e3
            sf_max = structureData[30:][30:].max()/Nballs/4
            ds = sf_max / resolution
            cax =ax.contourf(kx, ky, structureData/(Nballs), levels = np.arange(0,sf_max + ds,ds))
            ax.set_aspect('equal')
            cbar = fig.colorbar(cax)
            ax.set_title('Structure Factor density:'+str(density))
            ax.set_xlabel(r'$k_x (\frac{\pi}{a}$)')
            ax.set_ylabel(r'$k_y (\frac{\pi}{a}$)')

	
def plotDiffusionCoeff(rho,T):
    plt.figure()
    for temp in T:
        for density in rho:
            diffusionData = np.loadtxt('data2/ff-'+str(density)+'_temp-'+\
                str(temp)+'-dif.dat')
            #~ print diffusionData
            plt.semilogy(density,diffusionData,'k.')
            plt.title('Diffusion Coefficient v. Reduced Density')
            plt.xlabel(r'$\rho*$')
            plt.ylabel('D (length/iteration)')
            plt.grid()
            #~ plt.ylim(0)

parser = argparse.ArgumentParser(description='Which plots to make.')
parser.add_argument('-p','--plot', metavar='PLOT', action ="store",
                    default='diffusion', help='the plot to make')
parser.add_argument('-r','--rho', metavar = 'RHO',type = float, nargs = "+",
                    default = rhoDefault, help = "List of densities.")
parser.add_argument('-t','--temp', metavar = 'TEMP',type = float, nargs = "+",
                    default = TDefault, help = "List of Temperatures")
args = parser.parse_args()

print("\nPlot: %s"%args.plot)
print'Densities: ',args.rho
print'Temperatures: ',args.temp
		
if args.plot.lower() == 'diffusion':
        plotDiffusionCoeff(args.rho,args.temp)
elif args.plot.lower() == 'pressure':
    plotPressure(args.rho,args.temp)
elif args.plot.lower() == 'energy':
    energy(args.rho,args.temp)
elif (len(args.rho) <=5 and len(args.temp) <= 5) \
    and args.plot.lower()!=('diffusion' or 'pressure'):
    #~ if args.plot.lower() == 'energy':
        #~ plotEnergyPDF(args.rho,args.temp)
            
    if args.plot.lower() == 'radial':
        plotRadialDF(args.rho,args.temp)		

    elif args.plot.lower() == 'positions':
        plotPositions(args.rho,args.temp)
    elif args.plot.lower() == 'structure':
        plotStructureFactor(args.rho,args.temp)
elif (len(args.rho) > 5 or len(args.temp) > 5) \
    and args.plot.lower()!=('diffusion' or 'pressure'):
    print("\nPlease reduce the amount of plots you want to make by "\
        "shortening the density or temperature list to less than 5 each."\
        "\nYou can end up plotting way too many figures.\n")
    print("Number of plots asked for, Density: %d, Temp: %d, Total:"\
    "%d\n\n"%(len(args.rho),len(args.temp),len(args.rho)*len(args.temp)))
else:
    print("\nOr Plot type %s not recognized.\nPlease enter either:\n"\
                "diffusion\npositions\npressure\nor radial\n"%args.plot)

plt.show()
