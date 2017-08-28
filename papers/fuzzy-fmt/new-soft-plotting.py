#!/usr/bin/env python2
from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
########################################################################
## This script reads data from files and plots for new-soft
########################################################################
#Global Constants
x = 0
y = 1
z = 2

pmin = 1.9
pmax = 1.99
dp = 0.1

Tmin = 0.01
Tmax = 0.51
dT = 0.5

def plotPressure(pmin,pmax,dp,Tmin,Tmax,dT):
	nd = np.arange(pmin,pmax,dp)
	nT = np.arange(Tmin,Tmax,dT)
	pressure = np.zeros((len(nd),len(nT)))
	for density in nd:
		for temp in nT:
			pressure_file = open('data/ff-'+str(density)+'_temp-'+str(temp)+'-press.dat')
			with open('data/ff-'+str(density)+'_temp-'+str(temp)+'-press.dat') as f:
				for l in f:
					if(l[0] != "#"):
						tempPress = l.strip().split("\n")
						a = int(temp/dT)
						b = int((density/dp))-0

						pressure[b][a] = tempPress[0]
	fig = plt.figure()
	for i in range(len(nd)):
		print pressure[i][:]
		plt.plot(nT,pressure[i][:])
	plt.legend(nd[:])
	plt.title('Pressure v. Temperature')
	plt.xlabel('Temperature')
	plt.ylabel('Pressure')


def plotPositions(pmin,pmax,dp,Tmin,Tmax,dT):
	for density in np.arange(pmin,pmax,dp):
		for temp in np.arange(Tmin,Tmax,dT):
			position_file = open('data/ff-'+str(density)+'_temp-'+str(temp)+'-pos.dat','r')
			sphereNum = 0
			for line in position_file.readlines():
				if(line[0] != "#"):
					sphereNum += 1
			spheresx = np.zeros(sphereNum)
			spheresy = np.zeros(sphereNum)
			spheresz = np.zeros(sphereNum)
			count = 0
			with open('data/ff-'+str(density)+'_temp-'+str(temp)+'-pos.dat') as f:
				for l in f:
					if(l[0] != "#"):
						spherePos = l.strip().split("\t")
						spheresx[count] = spherePos[x]
						spheresy[count] = spherePos[y]
						spheresz[count] = spherePos[z]
						count += 1
			fig = plt.figure()
			ax = fig.add_subplot(111, projection = '3d')
			ax.scatter(spheresx,spheresy,spheresz,c = 'r', marker ='o')
			ax.set_title('Positions at Temp: '+ str(temp)+ ' and Density: ' + str(density))
			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')


def plotRadialDF(pmin,pmax,dp,Tmin,Tmax,dT):
	radboxes = np.zeros(1000)
	radheights = np.zeros(1000)
	for density in np.arange(pmin,pmax,dp):
		for temp in np.arange(Tmin,Tmax,dT):
			count = 0
			radial_file = open('data/ff-'+str(density)+'_temp-'+str(temp)+'-radial.dat')
			with open('data/ff-'+str(density)+'_temp-'+str(temp)+'-radial.dat') as f:
						for l in f:
							if(l[0] != "#"):
								data = l.strip().split("\t")
								radboxes[count] = data[0]
								radheights[count] = data[1]
								#~ radheights[count] = height[0]
			
			start = 300					count += 1
			f = np.fft.ifft(radheights[start:])

			plt.figure()
			plt.plot(radboxes[start:],radheights[start:])
			plt.title('Sum of Spheres at a Radial Distance, non-averaged. At temp: '+str(temp)+' and Density: ' +str(density))
			plt.xlabel('Radial Distance (r)')
			plt.ylabel('Number of Spheres at this distance')
			plt.figure()
			plt.plot(f)
			plt.title('FFT of Radial DF')

def plotEnergyPDF(pmin,pmax,dp,Tmin,Tmax,dT):
	for density in np.arange(pmin,pmax,dp):
		plt.figure()
		lis = []
		for temp in np.arange(Tmin,Tmax,dT):
			radial_file = open('data/ff-'+str(density)+'_temp-'+str(temp)+'-energy.dat')
			with open('data/ff-'+str(density)+'_temp-'+str(temp)+'-energy.dat') as f:
						for l in f:
							if(l[0] != "#"):
								data = l.strip().split("\t")
								mu = float(data[0])
								e2 = float(data[1])
			lis.append(str(density)+"-"+str(temp))
			sigma = np.sqrt(e2 - mu*mu)
			x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
			plt.plot(x,mlab.normpdf(x,mu,sigma))
			plt.title('Energy Probability Distribution Function Density: '+str(density))
			plt.xlabel('Total Internal Energy')
			plt.ylabel('Probability')
		plt.legend(labels = lis[:],title = 'Density-Temp',loc = 0)
	
def plotDiffusionCoeff(pmin,pmax,dp,Tmin,Tmax,dT):
	for density in np.arange(pmin,pmax,dp):
			plt.figure()
			for temp in np.arange(Tmin,Tmax,dT):	
				iterations = 0
				dif_file = open('data/ff-'+str(density)+'_temp-'+str(temp)+'-dif.dat')
				for line in dif_file.readlines():
					iterations += 1
				print "Iterations: ",iterations + 1
				D = np.zeros(iterations)
				count = 0
				with open('data/ff-'+str(density)+'_temp-'+str(temp)+'-dif.dat') as f:
							for l in f:
								if(l[0] != "#"):
									data = l.strip().split("\n")
									#~ print data
									D[count] = data[0]
									count += 1
				plt.plot(D)
				plt.title('Diffusion Coeff.')
				plt.xlabel('Iterations')
				plt.ylabel('Diffusion Coefficient')

#~ plotDiffusionCoeff(pmin,pmax,dp,Tmin,Tmax,dT)
#~ plotEnergyPDF(pmin,pmax,dp,Tmin,Tmax,dT)
plotRadialDF(pmin,pmax,dp,Tmin,Tmax,dT)
#~ plotPositions(pmin,pmax,dp,Tmin,Tmax,dT)		# Careful w/ this one
#~ plotPressure(pmin,pmax,dp,Tmin,Tmax,dT)

plt.show()
