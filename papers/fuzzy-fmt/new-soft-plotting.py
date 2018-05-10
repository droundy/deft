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
import new_soft_organize as nso
import csv
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit

#Global Constants
x = 0
y = 1
z = 2
direc = 'data/mc'

n_temps = nso.organize()[1]
temps = sorted(nso.organize()[0]['temp'])


def plot_diffusion(show,phase):

    '''Plots the Diffusion Coefficient v. Density.
    
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data. '''
    if phase:
        print 'Phase exists, showing Diffusion Coefficient'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Showing Diffusion Coefficient'
    plt.figure('''Diffusion Coefficient''')
    for t in temps:
        density = nso.grab_data(t)[0]
        diffusion = nso.grab_data(t)[1]

        plt.semilogy(density,diffusion,'.-',label = '%.1f'%(t))
        plt.legend(title = r'$T^*$')
        plt.title('Diffusion Coefficient v. Reduced Density')

        plt.xlabel(r'Reduced Density $(\rho^*)$')
        plt.ylabel(r'Diffusion Coefficient $\frac{r}{t_{max}}$')
        plt.ylim(0)

    if show == False:
        plt.savefig('figs/new-soft-diffusion.pdf')
        print 'Diffusion Coefficient Image Saved as figs/new-soft-diffusion.pdf'

def func(x,a,b,c):
    return a*x*x + b*x + c

def plot_pressure(show,phase):
    '''Plots the Pressure v. Density.
    
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data.'''
    
    if phase:
        print 'Phase exists, showing Pressure'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Showing Pressure'
        
    for t in temps:
        if t < 1.0:
            plt.figure('''Pressure v. Density Below 1.0''')
            density = nso.grab_data(t)[0]
            pressure = nso.grab_data(t)[2]
            plt.plot(density,pressure,'.-',label = '%.1f'%t)
            plt.legend(title = r'$T^*$')
            plt.title(r'Pressure v.Density $(T^* < 1.0)$')
            plt.ylabel(r'Reduced Pressure $(p^*)$')
            plt.xlabel(r'Reduced Density$(\rho^*)$')
            plt.grid()
            if show == False:
                plt.savefig('figs/new-soft-pressure_below.pdf')
                print 'Pressure 2 Image Saved as figs/new-soft-pressure_below.pdf'
        else:
            plt.figure('''Pressure v. Density Above 1.0''')
            density = nso.grab_data(t)[0]
            pressure = nso.grab_data(t)[2]
            plt.plot(density,pressure,'.-',label = '%.1f'%t)
            plt.legend(title = r'$T^*$')
            plt.title(r'Pressure v.Density $(T^* \geq 1.0)$')
            plt.ylabel(r'Reduced Pressure $(p^*)$')
            plt.xlabel(r'Reduced Density$(\rho^*)$')
            plt.grid()
            if show == False:
                plt.savefig('figs/new-soft-pressure_above.pdf')
                print 'Pressure 1 Image Saved as figs/new-soft-pressure_above.pdf'


def find_intersect(p_before,g_before,p_after,g_after):
    ''' Finds intersection of line segments for Gibbs Crossing Calculation'''
    for j in xrange(len(p_before)-1):
        p0 = np.array([p_before[j],g_before[j]])
        p1 = np.array([p_before[j+1],g_before[j+1]])
        for k in xrange(len(p_after)-1):
            q0 = np.array([p_after[k],g_after[k]])
            q1 = np.array([p_after[k+1],g_after[k+1]])
            params = np.linalg.solve(np.column_stack((p1-p0,q0-q1)),q0-p0)
            if (params[0] >=0) and (params[0]<=1) and (params[1] >=0 ) and (params[1] <= 1):
                return params[0]*(p1-p0)+p0
    
    
def free_energy(plot_volume,plot_helmholtz,plot_gibbs,show):
    ''' 1.  Solves for the Helmholtz Free Energy and Gibbs Free Energy
        2.  Determines the crossing point of the Gibbs to determine transitions
        3.  Plots a variety of interesting diagrams
        4.  Returns a dict containing pressures, densities, and volumes 
            corresponding to fluid and solid phase transitions '''
    phase = {'solid':{'T':[],'p':[],'d':[],'V':[],'dif':[]},'fluid':{'T':[],'p':[],'d':[],'V':[],'dif':[]}}
    for t in temps:
        density = nso.grab_data(t)[0]
        pressure = nso.grab_data(t)[2]
        volume = nso.grab_data(t)[3]
        Nballs = nso.grab_data(t)[4]
        diffusion = nso.grab_data(t)[1]

        if plot_volume:
            plt.figure('''Pressure v. Volume''')
            plt.plot(volume,pressure, '.-', label = '%.1f'%t)
            plt.legend(title = r'$T^*$')
            plt.title('Pressure v.Volume')
            plt.ylabel('Pressure')
            plt.xlabel('Volume')
            plt.grid()
            if show == False:
                plt.savefig('figs/new-soft-pressure_v_volume.pdf')
                print 'Pressure v Volume Image Saved as figs/new-soft-pressure_v_volume.pdf'
        
        hfree_length = len(volume)
        hfree = np.zeros(hfree_length)

        for i in range(1,hfree_length):
            hfree[i] = hfree[i-1] + (pressure[i]+pressure[i-1])*(volume[i]-volume[i-1])/(2*Nballs)
        hfree = -hfree
        ''' This section is curve fitting the last temperature to add to my thesis '''
        lower_curve = np.zeros(hfree_length)
        upper_curve = np.zeros(hfree_length)

        if t == 10.0:
            for j in range(1,hfree_length):
                if pressure[j-1] > pressure[j]:
                    break
            lower_opt, lower_cov = curve_fit(func,density[:j],hfree[:j])
            upper_opt, upper_cov = curve_fit(func,density[j:],hfree[j:])
            
            for j in range(hfree_length):
                lower_curve[j] = func(density[j], lower_opt[0], lower_opt[1], lower_opt[2])
                upper_curve[j] = func(density[j], upper_opt[0], upper_opt[1], upper_opt[2])

            ## Working on the Curve fitting for the free energy
        if plot_helmholtz:
            plt.figure('''Helmholtz Free Energy''')
            plt.plot(density,hfree,'.-', label = '%.1f'%t, markersize = 2)
            plt.title('Reduced Helmholtz Free Energy v. Density',fontsize = "x-large")
            plt.xlabel(r'Reduced Density $(\rho^*)$',fontsize = "x-large")
            plt.ylabel(r'Reduced Helmholtz Free Energy $(f - f_0)$',fontsize = "x-large")
            plt.grid()
            #~ if t == 10.0:
                #~ plt.plot(density[:j], lower_curve[:j], color = 'r', ls = '--', lw = 0.5, label = 'Fluid Line')
                #~ plt.plot(density[:j], upper_curve[:j], color = 'g', ls = '--', lw = 0.5, label = 'Solid Line')
            plt.legend(title = r'$T^*$')
            if show == False:
                plt.savefig('figs/new-soft-helmholtz.pdf')
                print 'Helmholtz Free Energy Image Saved as figs/new-soft-helmholtz.pdf'
            
        gibbs = np.zeros(hfree_length)

        for i in range(hfree_length):
            gibbs[i] = hfree[i] + pressure[i]*volume[i]/Nballs

        
        direction_change = [-1,-1]
        ''' Solves for Crossing point of Gibbs v. Pressure '''
        for i in range(len(gibbs)-1):
            greater = gibbs[i+1] > gibbs[i]
            less= gibbs[i+1] < gibbs[i]
            if (less) and (direction_change[0] == -1):
                direction_change[0] = i + 1
            elif (greater) and (direction_change[0] != -1) :
                direction_change[1] = i + 1
                break
        if direction_change[0] != -1:
            p_before = pressure[:direction_change[0]]
            g_before = gibbs[:direction_change[0]]
            p_after = pressure[direction_change[1]:]
            g_after = gibbs[direction_change[1]:]
            cross = find_intersect(p_before,g_before,p_after,g_after)
            
            
            gdiff = [max(gibbs)-min(gibbs),max(gibbs)-min(gibbs)]
            nearest = [-1,-1] 
            
            for i in range(len(g_before)):
                tgdiff = g_before[i]-cross[1]
                if abs(tgdiff) < abs(gdiff[0]):
                    nearest[0],gdiff[0] = i,tgdiff

            for j in range(len(g_after)):
                tgdiff = g_after[j]-cross[1]
                if abs(tgdiff) < abs(gdiff[1]):
                    nearest[1],gdiff[1] = (len(gibbs)-len(g_after)+j),tgdiff

            phase['solid']['T'].append(t)
            phase['solid']['p'].append(pressure[max(nearest)])
            phase['solid']['V'].append(volume[max(nearest)])
            phase['solid']['d'].append(density[max(nearest)])
            phase['solid']['dif'].append(diffusion[max(nearest)])
            
            phase['fluid']['T'].append(t)
            phase['fluid']['p'].append(pressure[min(nearest)])
            phase['fluid']['V'].append(volume[min(nearest)])
            phase['fluid']['d'].append(density[min(nearest)])
            phase['fluid']['dif'].append(diffusion[min(nearest)])
            
            

            
            if plot_gibbs:
                plt.figure()
                plt.axhline(y = cross[1], color = '0.3', ls = ':', lw = '0.5')
                plt.axvline(x = cross[0], color = '0.3', ls = ':', lw = '0.5')
                plt.plot(pressure,gibbs,'k.-')
                plt.plot(cross[0],cross[1],'r.')
                plt.plot(pressure[min(nearest)],gibbs[(min(nearest))], 'm.')
                plt.plot(pressure[max(nearest)],gibbs[(max(nearest))], 'm.')
                #~ plt.annotate(s ='Cross at (p=%.2f,g=%.2f), Nearest Densities: %s, %s'
                    #~ %(cross[0],cross[1],density[nearest[0]],density[nearest[1]]),
                    #~ xy = [cross[0]+0.1*t,cross[1]-0.1*t])
                plt.title(r'Gibbs Free Energy $(T^*$ = %.1f)'%t,fontsize = "x-large")
                plt.xlabel(r'Reduced Pressure $(p^*)$',fontsize = "x-large")
                plt.ylabel(r'Reduced Gibbs Free Energy $(g - g_0)$)',fontsize = "x-large")
                plt.grid()
                if show == False:
                    plt.savefig('figs/new-soft-gibbs_T-%s.pdf'%t)
                    print 'Gibbs Free Energy Image Saved as figs/new-soft-gibbs_T-%s.pdf'%t
                else:
                    print("Gibbs Crossing Point at T: %s, Pressure: %s, Gibbs: %s")\
                        %(t,cross[0],cross[1])
    return phase



def plot_PT_diagram(show,phase):
    ''' Plots the Pressure-Temperature Phase Diagram.
    
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data.'''
    if phase:
        print 'Phase exists, showing Pressure-Temperature Phase Diagram'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Showing Pressure-Temperature Phase Diagram'
        phase = free_energy(False,False,False,False)
        
    plt.figure('''PT Diagram''')
    plt.plot(phase['fluid']['T'],phase['fluid']['p'],'r.-', label = 'Fluid')
    plt.fill_between(phase['fluid']['T'],0,phase['fluid']['p'], color = 'red', alpha = '0.5')
    plt.plot(phase['solid']['T'],phase['solid']['p'],'b.-', label = 'Solid')
    plt.fill_between(phase['solid']['T'],phase['solid']['p'],max(phase['solid']['p']), color = 'blue', alpha = '0.5')
    plt.title('Freezing and Melting Points  (P-T Diagram)')
    plt.xlabel(r'Reduced Temperature $(T^*)$')
    plt.ylabel(r'Reduced Pressure $(p^*)$')
    plt.legend()

    if show == False:
        plt.savefig('figs/new-soft-PT_diagram.pdf')
        print 'pT Phase Diagram Image Saved as figs/new-soft-PT_diagram.pdf'

def plot_prho_diagram(show,phase):
    ''' Plots the Pressure-Density Phase Diagram.
    
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data.'''
    if phase:
        print 'Phase exists, showing Pressure-Density Phase Diagram'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Showing Pressure-Density Phase Diagram'
        phase = free_energy(False,False,False,False)
        
    plt.figure('''p-rho Diagram''')
    plt.plot(phase['fluid']['d'],phase['fluid']['p'],'r.-', label = 'Fluid')
    plt.fill_between(phase['fluid']['d'],phase['fluid']['p'],max(phase['fluid']['p']), color = 'red', alpha = '0.5')
    plt.plot(phase['solid']['d'],phase['solid']['p'],'b.-', label = 'Solid')
    plt.fill_between(phase['solid']['d'],0,phase['solid']['p'], color = 'blue', alpha = '0.5')
    plt.title(r'Freezing and Melting Points  $(p^*-\rho$ Diagram)')
    plt.xlabel(r'Reduced Density $(\rho^*)$')
    plt.ylabel(r'Reduced Pressure $(p^*)$')
    plt.legend()
    if show == False:
        plt.savefig('figs/new-soft-Trho_diagram.pdf')
        print 'TRho Phase Diagram Image Saved as figs/new-soft-Trho_diagram.pdf'


def phase_table(phase):
    ''' Prints out a table of densities and pressures corresponding to 
        phase transitions 
        
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data.
            '''
    if phase:
        print 'Phase exists, showing Data Table'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Showing Data Table'
        phase = free_energy(False,False,False,False)

    print 'Solid Temperatures: ',phase['solid']['T']
    print 'Solid Pressures: ',phase['solid']['p']
    print 'Solid Densities: ',phase['solid']['d']
    print 'Solid Diffusion: ',phase['solid']['dif']


    print 'Fluid Temperatures: ',phase['fluid']['T']
    print 'Fluid Pressures: ',phase['fluid']['p']
    print 'Fluid Densities: ',phase['fluid']['d']
    print 'Fluid Diffusion: ',phase['fluid']['dif']


def plotRadialDF(n,show,phase):
    ''' Plots the fluid and solid radial distribution function for the 
        fluid and solid as determined by the free_energy function.
        
        "n" is the index of the temperature in phase that will be plotted
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data. '''
    if phase:
        print 'Phase exists, plotting Radial Distribution Function.'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Plotting Radial Distribution Function'
        phase = free_energy(False,False,False,False)
    TSolid = phase['solid']['T'][n]
    TFluid = phase['fluid']['T'][n]
    dSolid = phase['solid']['d'][n]
    dFluid = phase['fluid']['d'][n]
    radialDataSolid = np.loadtxt('%s/ff-%s_temp-%s-%s'%(direc,dSolid,TSolid,'radial.dat'))
    radialDataFluid = np.loadtxt('%s/ff-%s_temp-%s-%s'%(direc,dFluid,TFluid,'radial.dat'))
    plt.figure('''Radial Distribution Function ''')
    plt.plot(radialDataSolid[:,0],2*radialDataSolid[:,1]/256,'k', label = r'Solid: $\rho$=%.2f'%dSolid)
    plt.plot(radialDataFluid[:,0],2*radialDataFluid[:,1]/256,'k--', label = r'Fluid: $\rho$=%.2f'%dFluid)
    plt.axvline(x = 2**(1/6), color = '0.5', ls = '-.', label = 'WCA Minimum')
    plt.title(r'Radial Distribution Function $(T^*$ =%.2f )'%TSolid)
    plt.xlabel(r'Radial Distance $(\frac{r}{\sigma})$')
    plt.ylabel(r'$g(r)$')
    plt.legend(title = r'Reduced Density $(\rho^*)$')
    plt.grid()
    if show == False:
        plt.savefig('figs/new-soft-rdf_T-%.2f.pdf'%TFluid)
        print 'Radial Distribution Function Image Saved as figs/new-soft-rdf_T-%.2f.pdf'%TFluid


def plotPositions(rho,T):
    ''' Plots a 3D image to show spheres in real space '''
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


def plotStructureFactor(n,show,phase):
    '''Plots the structure factor.
    
        "n" is the index of the temperature in phase that will be plotted
        "show" determines whether or not to save the image to file
        "phase" is where all the file information is stored,
            a value of False in here will run the free_energy() function
            to obtain the data.'''
    if phase:
        print 'Phase exists, plotting Structure Factor. This takes the longest.'
    else:
        print 'Phase does not exist, running free_energy()'
        print '\t Plotting Structure Factor'
        phase = free_energy(False,False,False,False)
        
    TSolid = phase['solid']['T'][n]
    TFluid = phase['fluid']['T'][n]
    dSolid = phase['solid']['d'][n]
    dFluid = phase['fluid']['d'][n]
    Nballs = nso.grab_data(0.5)[4]
    fnameFluid = '%s/ff-%s_temp-%s-%s'%(direc,dFluid,TFluid,'struc.dat')
    fnameSolid = '%s/ff-%s_temp-%s-%s'%(direc,dSolid,TSolid,'struc.dat')

    structureDataFluid = np.loadtxt(fnameFluid)
    structureDataSolid = np.loadtxt(fnameSolid)
    dk = 1/40

    kvalsSolid = dk*np.arange(0.0, len(structureDataSolid), 1)
    kxSolid, kySolid = np.meshgrid(kvalsSolid, kvalsSolid)

    kvalsFluid = dk*np.arange(0.0, len(structureDataFluid), 1)
    kxFluid, kyFluid = np.meshgrid(kvalsFluid, kvalsFluid)

    resolution = 2e3
    sf_max_Fluid = structureDataFluid[30:][30:].max()/Nballs/2
    sf_max_Solid = structureDataSolid[30:][30:].max()/Nballs/2

    ds_Fluid = sf_max_Fluid / resolution
    ds_Solid = sf_max_Solid / resolution

    fig = plt.figure('''Structure Factor 1''')
    #~ st = fig.suptitle(r'Fluid and Solid Structure Factors $(T^*$ = %.2f)' %(TFluid), fontsize = "x-large")
    #~ plt.subplot(2,1,1, aspect = 'equal')
    cset1 = plt.contourf(kxFluid, kyFluid, structureDataFluid/(Nballs), 
        levels = np.arange(0,sf_max_Fluid + ds_Fluid,ds_Fluid), aspect = 'equal')
    #~ plt.colorbar(cset1)
    plt.title('Fluid at Density: %.2f'%(dFluid))
    plt.xlabel(r'$k_x (\frac{\pi}{a}$)', fontsize = "x-large")
    plt.ylabel(r'$k_y (\frac{\pi}{a}$)', fontsize = "x-large")
    plt.tight_layout()
    if show == False:
        plt.savefig('figs/new-soft-sf_T-%.2f_d-%.2f.png'%(TFluid,dFluid))
        print 'Structure Factor Image Saved as figs/new-soft-sf_T-%.2f.png'%TFluid
    
    fig2 = plt.figure('''Structure Factor 2''')
    #~ plt.subplot(2,2,1, aspect = 'equal')
    cset2 = plt.contourf(kxSolid, kySolid, structureDataSolid/(Nballs), 
        levels = np.arange(0,sf_max_Solid + ds_Solid,ds_Solid),aspect = 'equal')
    #~ plt.colorbar(cset2)
    plt.title('Solid at Density: %.2f'%(dSolid))
    plt.xlabel(r'$k_x (\frac{\pi}{a}$)', fontsize = "x-large")
    plt.ylabel(r'$k_y (\frac{\pi}{a}$)', fontsize = "x-large")
    if show == False:
        plt.savefig('figs/new-soft-sf_T-%.2f_d-%.2f.png'%(TFluid,dSolid))
        print 'Structure Factor Image Saved as figs/new-soft-sf_T-%.2f.png'%TFluid
        # Cannot be saved as a PDF, too much stuff.

print("\nPlot: %s"%args.plot)
		
if args.plot.lower() == '':
    # n = 3 for 1.0, n = 7 for 10.0
    n = 7
    show = True
    plot_volume=False
    plot_helmholtz=False
    plot_gibbs=True
    phase = free_energy(plot_volume,plot_helmholtz,plot_gibbs,show)
    #~ phase_table(phase)
    plotRadialDF(n,show,phase)
    #~ plotStructureFactor(n,show,phase)
    #~ plot_diffusion(show,phase)
    #~ plot_pressure(show,phase)
    #~ plot_PT_diagram(show,phase)
    #~ plot_prho_diagram(show,phase)
    
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


if show:
    plt.show()
