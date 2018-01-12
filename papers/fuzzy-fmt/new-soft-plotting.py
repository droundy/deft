#!/usr/bin/env python2
from __future__ import division
import argparse

parser = argparse.ArgumentParser(description='Which plots to make.')
parser.add_argument('-p','--plot', metavar='PLOT', action ="store",
                    default='', help='the plot to make')
parser.add_argument('--no-show', help='do not pop up windows',
                    action="store_true")
parser.add_argument('-r','--rho', metavar = 'RHO',type = float, nargs = "+",
                    default = [1.0], help = "List of densities.")
parser.add_argument('-t','--temp', metavar = 'TEMP',type = float, nargs = "+",
                    default = [1.0], help = "List of Temperatures")
args = parser.parse_args()

if args.no_show:
    # We need the following two lines in order for matplotlib to work
    # without access to an X server.
    import matplotlib
    matplotlib.use('Agg')

from mpl_toolkits.mplot3d import Axes3D
import csv
import numpy as np
import matplotlib.pyplot as plt
import glob

#Global Constants
x = 0
y = 1
z = 2

direc = 'data/mc'
cpalette = ['b','g','r','c','m','y','k']

filenames = glob.glob('%s/*.dat'%(direc))
file_dict = {'temp':{}}
temps,ptypes,dens = [],[],[]

for f in filenames:
	''' Finds all the densities, temperatures, and plot types from files'''
	t,p,d = f.split('-')[2],f.split('-')[3],f.split('-')[1].split('_')[0]
	temps.append(t),ptypes.append(p),dens.append(d)
	temps,ptypes = np.unique(temps).tolist(),np.unique(ptypes).tolist()
	dens = np.unique(dens).tolist()
for t in temps:
	'''Creates a dictionary of all available data'''
	file_dict['temp'].update({t:{}})
	for p in ptypes:
		file_dict['temp'][t].update({p:[]})
		for d in dens:
			if ('%s/ff-%s_temp-%s-%s'%(direc,d,t,p)) in filenames:
				file_dict['temp'][t][p].append(d)

def all_data():
    ''' Plots many Diffusion and Pressure '''
    density,diffusion_data,pressure_data,volume = {},{},{},{}
    nsubplot = len(file_dict['temp'])
    for t in sorted(file_dict['temp']):
        color = file_dict['temp'].keys().index(t)
        density.update({t:[]})
        diffusion_data.update({t:[]})
        pressure_data.update({t:[]})
        volume.update({t:[]})
        for d in file_dict['temp'][t]['dif.dat']:
            density[t].append(float(d))
            diffusion_data[t].append(np.loadtxt('%s/ff-%s_temp-%s-%s'%
                (direc,d,t,'dif.dat')))
            pressure_data[t].append(np.loadtxt('%s/ff-%s_temp-%s-%s'%
                (direc,d,t,'press.dat')))
            with open('%s/ff-%s_temp-%s-%s'%(direc,d,t,'press.dat')) as f:
                for l in f.readlines():
                    if l[:10] == '# volume: ':
                        volume[t].append(float(l[10:]))
                        break
            Nballs = -1
            with open('%s/ff-%s_temp-%s-%s'%(direc,d,t,'press.dat')) as f:
                for l in f.readlines():
                    if l[:5] == '# N: ':
                        NBalls = int(l[5:])
                        break
        F = np.zeros(len(volume[t])-1)
        p,v = np.zeros(len(F)),np.zeros(len(F))
        F[0] = -(pressure_data[t][1]+pressure_data[t][0])*\
            (volume[t][1]-volume[t][0])/(2*Nballs)
        v[0] = (volume[t][1]+volume[t][0])/(2*Nballs)
        p[0] = (pressure_data[t][1]+pressure_data[t][0])/2
        for n in range(1,len(F)):
            F[n] += F[n-1]-(pressure_data[t][n+1]+pressure_data[t][n])*\
                (volume[t][n+1]-volume[t][n])/(2*Nballs)
            v[n] = (volume[t][n+1]+volume[t][n])/(2*Nballs)
            p[n] = (pressure_data[t][n+1]+pressure_data[t][n])/2
        gibbs = F + p*v 

        direction_change = [-1,-1]
        ''' Solves for Crossing point of Gibbs v. Pressure '''
        for i in range(len(gibbs)-1):
            greater = gibbs[i+1] > gibbs[i]
            less= gibbs[i+1] < gibbs[i]
            if (greater) and (direction_change[0] == -1):
                direction_change[0] = i
            elif (less) and (direction_change[0] != -1) :
                direction_change[1] = i
                break
        p_before = p[:direction_change[0]]
        g_before = gibbs[:direction_change[0]]
        p_after = p[direction_change[1]:]
        g_after = gibbs[direction_change[1]:]
        cross = find_intersect(p_before,g_before,p_after,g_after)
        
        gdiff = [max(gibbs)-min(gibbs),max(gibbs)-min(gibbs)]
        nearest = [-1,-1] # Two closest indices to the intersection

        ''' Diffusion Coefficient'''
        #~ plt.figure(1)
        #~ plt.semilogy(density[t],diffusion_data[t],'%s.-'%(cpalette[color]))
        #~ plt.legend(sorted(file_dict['temp']))
        #~ plt.title('Diffusion Coefficient v. Reduced Density')
        #~ plt.xlabel(r'$\rho*$')
        #~ plt.ylabel('D (length/iteration)')
        #~ plt.ylim(0)
        #~ ''' Pressure v. Density'''
        #~ plt.figure(2)
        #~ plt.plot(density[t],pressure_data[t],'%s.-'%(cpalette[color]))
        #~ plt.legend(sorted(file_dict['temp']))
        #~ plt.title('Pressure v.Density')
        #~ plt.ylabel('Pressure')
        #~ plt.xlabel('Density')
        #~ plt.grid()
        #plt.savefig('figs/FIXME.pdf')
        #~ ''' Pressure v. Volume '''
        #~ plt.figure(3)
        #~ plt.plot(volume[t],pressure_data[t],'%s.-'%(cpalette[color]))
        #~ plt.xlabel('Volume')
        #~ plt.ylabel('Pressure')
        #~ plt.title('Pressure v. Volume')
        #~ plt.legend(sorted(file_dict['temp']))
        #~ plt.grid()
        #~ ''' Helmholtz Free Energy'''
        #~ plt.figure(4)
        #~ plt.plot(v,F,'%s.-'%(cpalette[color]))
        #~ plt.title(r'$f - f_0$')
        #~ plt.xlabel('reduced volume v')
        #~ plt.ylabel('reduced helmholtz')
        #~ plt.legend(sorted(file_dict['temp']))
        #~ plt.grid()

        if cross != None:
            for i in range(len(g_before)):
                tgdiff = g_before[i]-cross[1]
                if abs(tgdiff) < gdiff[0]:
                    nearest[0],gdiff[0] = i,tgdiff

            for j in range(len(g_after)):
                tgdiff = g_after[j]-cross[1]
                if abs(tgdiff) < gdiff[1]:
                    nearest[1],gdiff[1] = (len(gibbs)-len(g_after)+j),tgdiff

            #~ ''' Phase Diagram '''
            plt.figure(5)
            #~ print min(nearest),max(nearest)
            plt.plot(density[t][min(nearest)],float(t),'r.')
            plt.plot(density[t][max(nearest)],float(t),'b.')
            plt.title('Freezing and Melting Points')
            plt.xlabel('Reduced Density')
            plt.ylabel('Reduced Temperature')
            ''' Gibbs Free Energy '''
            plt.figure()
            plt.plot(p,gibbs,'k.-')
            plt.plot(cross[0],cross[1],'r.')
            plt.plot(p[nearest],gibbs[nearest],'m.')
            plt.annotate(s ='Cross at (p=%.2f,g=%.2f), Nearest Densities: %s, %s'
                %(cross[0],cross[1],density[t][nearest[0]],density[t][nearest[1]]),
                xy = [cross[0],cross[1]])
            plt.xlabel('pressure p*')
            plt.ylabel('g(v)')
            plt.title(r' Reduced Gibbs Temp: %s'%t)
            plt.grid()
            ''' Helmholtz Free w/ Chemical Potential '''
            #~ plt.figure()
            #~ mu = -72
            #~ rho = np.array(pressure_data[t][:-1])
            #~ plt.plot(rho,F - mu*rho)
            #~ plt.xlabel('Density')
            #~ plt.ylabel(r'F - $\mu*\rho$')
            #~ plt.title('Helmholtz Energy and Chemical Potential')

    
def find_intersect(p_before,g_before,p_after,g_after):
    ''' Finds intersection of line segments '''
    for j in xrange(len(p_before)-1):
        p0 = np.array([p_before[j],g_before[j]])
        p1 = np.array([p_before[j+1],g_before[j+1]])
        for k in xrange(len(p_after)-1):
            q0 = np.array([p_after[k],g_after[k]])
            q1 = np.array([p_after[k+1],g_after[k+1]])
            params = np.linalg.solve(np.column_stack((p1-p0,q0-q1)),q0-p0)
            if (params[0] >=0) and (params[0]<=1) and (params[1] >=0 ) and (params[1] <= 1):
                return params[0]*(p1-p0)+p0
            #~ else:
                #~ return [-1,-1]

        

def plotPositions(rho,T):
	for d in rho:
		for t in T:
			spheres = np.loadtxt('%s/ff-%s_temp-%s-%s'%(direc,d,t,'pos.dat'))
			fig = plt.figure()
			ax = fig.add_subplot(111, projection = '3d')
			ax.scatter(spheres[:,0],spheres[:,1],spheres[:,2],\
                c = 'r', marker ='o')
			ax.set_title('Positions at Temp: '+ str(t)+\
                ' and Density: ' + str(d))
			ax.set_xlabel('x')
			ax.set_ylabel('y')
			ax.set_zlabel('z')


def plotRadialDF(rho,T):
	for d in rho:
		for t in T:
			radialData  = np.loadtxt('%s/ff-%s_temp-%s-%s'%(direc,d,t,'radial.dat'))
			plt.figure()
			plt.plot(radialData[:,0],2*radialData[:,1]/256)#/(radialData[:,0]**2))
			plt.title('Radial Distribution Function.\n'\
				'Temp: '+str(t)+' and Density: ' +str(d))
			plt.xlabel('Radial Distance (r)')
			plt.ylabel('Nearby Spheres')
                plt.grid()

        
def plotStructureFactor(rho,T):
    for t in T:
        for d in rho:
            fig, ax = plt.subplots()
            fname = '%s/ff-%s_temp-%s-%s'%(direc,d,t,'struc.dat')
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
            sf_max = structureData[30:][30:].max()/Nballs/2
            ds = sf_max / resolution
            cax =ax.contourf(kx, ky, structureData/(Nballs), levels = np.arange(0,sf_max + ds,ds))
            ax.set_aspect('equal')
            cbar = fig.colorbar(cax)
            ax.set_title('Structure Factor density:'+str(d))
            ax.set_xlabel(r'$k_x (\frac{\pi}{a}$)')
            ax.set_ylabel(r'$k_y (\frac{\pi}{a}$)')


print("\nPlot: %s"%args.plot)
		
if args.plot.lower() == '':
        all_data()
elif args.plot.lower() == 'energy':
    energy(args.rho,args.temp)
elif (len(args.rho) <=5 and len(args.temp) <= 5) \
    and args.plot.lower()!=(''):
    print'Densities: ',args.rho
    print'Temperatures: ',args.temp
    if args.plot.lower() == 'radial':
        plotRadialDF(args.rho,args.temp)		

    elif args.plot.lower() == 'positions':
        plotPositions(args.rho,args.temp)
    elif args.plot.lower() == 'structure':
        plotStructureFactor(args.rho,args.temp)
elif (len(args.rho) > 5 or len(args.temp) > 5) \
    and args.plot.lower()!=(''):
    print("\nPlease reduce the amount of plots you want to make by "\
        "shortening the density or temperature list to less than 5 each."\
        "\nYou can end up plotting way too many figures.\n")
    print("Number of plots asked for, Density: %d, Temp: %d, Total:"\
    "%d\n\n"%(len(args.rho),len(args.temp),len(args.rho)*len(args.temp)))
else:
    print("\nOr Plot type %s not recognized.\nPlease enter either:\n"\
                "positions\nstructure\nor radial\n"%args.plot)

plt.show()
