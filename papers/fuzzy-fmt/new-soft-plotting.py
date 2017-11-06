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
rhoDefault = [0.7]
TDefault = [2.0]

def plotPressure(rho,T):
	plt.figure()
	for density in rho:
		pressArray = []
		for temp in T:
			pressure = np.loadtxt('data/ff-'+str(density)+\
                '_temp-'+str(temp)+'-press.dat')
                
                    
                pressArray.append(float(pressure))
                plt.plot(density,pressArray,'k.')
        plt.title('Pressure-Density at T*: '+str(temp))
        plt.ylabel('Pressure p*')
        plt.xlabel(r'Density $\rho*$')

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


#~ def plotEnergyPDF(rho,T):
	#~ for density in rho:
		#~ plt.figure()
		#~ for temp in T:
			#~ energyData = np.loadtxt('data/ff-'+str(density)+'_temp-'+\
                #~ str(temp)+'-energy.dat')
			#~ mu,e2 = energyData
			#~ sigma = np.sqrt(e2 - mu*mu)
			#~ x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
			#~ plt.plot(x,mlab.normpdf(x,mu,sigma))
			#~ plt.title('Energy Probability Distribution Function Density: '+\
                #~ str(density))
			#~ plt.xlabel('Total Internal Energy')
			#~ plt.ylabel('Probability')
            
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
            diffusionData = np.loadtxt('data/ff-'+str(density)+'_temp-'+\
                str(temp)+'-dif.dat')
            print diffusionData
            plt.semilogy(density,diffusionData,'k.')
            plt.title('Diffusion Coefficient v. Reduced Density')
            plt.xlabel(r'$\rho*$')
            plt.ylabel('D (length/iteration)')
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
