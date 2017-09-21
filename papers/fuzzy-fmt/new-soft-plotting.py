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
# Will crash if Position, Energy, or Radial are called.
rhoDefault = [0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,1.0,1.1,1.2,1.3,1.4]
TDefault = [1.0]

def plotPressure(rho,T):
	plt.figure()
	for density in rho:
		pressArray = []
		for temp in T:
			pressure = np.loadtxt('data/ff-'+str(density)+\
                '_temp-'+str(temp)+'-press.dat')
			pressArray.append(float(pressure))
                plt.plot(pressArray,density,'k.-')
        plt.title('Density v. Pressure')
        plt.xlabel('Pressure')
        plt.ylabel('Density')

def plotPositions(rho,T):
	for density in rho:
		for temp in T:
			spheres = np.loadtxt('data/ff-'+str(density)+\
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
			radialData  = np.loadtxt('data/ff-'+str(density)+\
                '_temp-'+str(temp)+'-radial.dat')
			plt.figure()
			plt.plot(radialData[:,0],radialData[:,1]/(radialData[:,0]**2))
			plt.title('Sum of Spheres at a Radial Distance, non-averaged.\n'\
				'Temp: '+str(temp)+' and Density: ' +str(density))
			plt.xlabel('Radial Distance (r)')
			plt.ylabel('Number of Spheres at this distance')


def plotEnergyPDF(rho,T):
	for density in rho:
		plt.figure()
		for temp in T:
			energyData = np.loadtxt('data/ff-'+str(density)+'_temp-'+\
                str(temp)+'-energy.dat')
			mu,e2 = energyData
			sigma = np.sqrt(e2 - mu*mu)
			x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
			plt.plot(x,mlab.normpdf(x,mu,sigma))
			plt.title('Energy Probability Distribution Function Density: '+\
                str(density))
			plt.xlabel('Total Internal Energy')
			plt.ylabel('Probability')

	
def plotDiffusionCoeff(rho,T):
    plt.figure()
    for temp in T:
        for density in rho:
            diffusionData = np.loadtxt('data/ff-'+str(density)+'_temp-'+\
                str(temp)+'-dif.dat')
            plt.semilogy(density,diffusionData[999,1],'ko')
            plt.title('Diffusion Coefficient v. Reduced Density')
            plt.xlabel(r'$\rho*$')
            plt.ylabel('D (length/iteration)')

parser = argparse.ArgumentParser(description='Which plots to make.')
parser.add_argument('--plot', metavar='PLOT', action ="store",
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
    if args.plot.lower() == 'energy':
        plotEnergyPDF(args.rho,args.temp)
            
    elif args.plot.lower() == 'radial':
        plotRadialDF(args.rho,args.temp)		

    elif args.plot.lower() == 'positions':
        plotPositions(args.rho,args.temp)
elif (len(args.rho) > 5 or len(args.temp) > 5) \
    and args.plot.lower()!=('diffusion' or 'pressure'):
    print("\nPlease reduce the amount of plots you want to make by "\
        "shortening the density or temperature list to less than 5 each."\
        "\nYou can end up plotting way too many figures.\n")
    print("Number of plots asked for, Density: %d, Temp: %d, Total:"\
    "%d\n\n"%(len(args.rho),len(args.temp),len(args.rho)*len(args.temp)))
else:
    print("\nOr Plot type %s not recognized.\nPlease enter either:\n"\
                "diffusion\nenergy\npositions\npressure\nor radial\n"%args.plot)

plt.show()
