from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import csv
import numpy as np
import matplotlib.pyplot as plt

x = 0
y = 1
z = 2
arraylength = 0

ffmin = 0.1
ffmax = 1.1
dff = 0.1

tempmin = 0.01
tempmax = 3.5
dtemp = 0.5


tempFF = np.arange(ffmin,ffmax,dff)
tempTemp = np.arange(tempmin,tempmax,dtemp)

pressure = np.zeros((len(tempFF),len(tempTemp)))


#######################################################################
## Make sure that the min and max for the ff and temperature are the same
## as in the Chris-calls.py file or else they won't match :)
#######################################################################

for fillFrac in tempFF:
	print fillFrac
	for temp in tempTemp:
		print temp
		position_file = open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-pos.dat')
		radial_file = open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-radial.dat')
		pressure_file = open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-press.dat')
		
		sphereNum = 0
		for line in position_file.readlines():
			sphereNum += 1

		print "You could try...", np.loadtxt('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-press.dat')
		with open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-press.dat') as f:
			for l in f:
				if(l[0] != "#"):
					tempPress = l.strip().split("\n")
					a = int(temp/dtemp)
					b = int(fillFrac/dff)-1
					pressure[b][a] = tempPress[0]
				
		
		#~ spheresx = np.zeros(sphereNum)
		#~ spheresy = np.zeros(sphereNum)
		#~ spheresz = np.zeros(sphereNum)
		
		#~ radboxes = np.linspace(0,2**(5/6),1000)
		#~ radheights = np.zeros(1000)
		#~ print temp
		
		#~ count = 0
		#~ with open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-pos.dat') as f:
			#~ for l in f:
				#~ if(l[0] != "#"):
					#~ spherePos = l.strip().split("\t")
					#~ spheresx[count] = spherePos[x]
					#~ spheresy[count] = spherePos[y]
					#~ spheresz[count] = spherePos[z]
					#~ count += 1
		
		#~ count = 0
		#~ with open('CHRIS/ff-'+str(fillFrac)+'_temp-'+str(temp)+'-radial.dat') as f:
			#~ for l in f:
				#~ if(l[0] != "#"):
					#~ height = l.strip().split("\n")
					#~ radheights[count] = height[0]
					#~ count += 1
					
		#~ fig = plt.figure()
		#~ ax = fig.add_subplot(111, projection = '3d')
		#~ ax.scatter(spheresx,spheresy,spheresz,c = 'r', marker ='o')
		#~ ax.set_title('Positions at Temp: '+ str(temp)+ ' and FF: ' + str(fillFrac))
		#~ ax.set_xlabel('x')
		#~ ax.set_ylabel('y')
		#~ ax.set_zlabel('z')
		
		#~ plt.figure()
		#~ plt.plot(radboxes,radheights)
		#~ plt.title('Sum of Spheres at a Radial Distance, non-averaged. At temp: '+str(temp)+' and FF: ' +str(fillFrac))
		#~ plt.xlabel('Radial Distance (r)')
		#~ plt.ylabel('Number of Spheres at this distance')
	
print pressure

plt.figure()
for i in range(len(tempFF)):
	plt.plot(tempTemp,pressure[i][:])
plt.title('Pressure v. Temperature')
plt.xlabel('Temperature')
plt.ylabel('Pressure')


plt.show()
