#!/usr/bin/python3
from __future__ import division
import matplotlib as mp
#mp.use('Agg')
from pylab import *
import os, sys, glob, socket, argparse

import args_figs

# parse arguments
parser = argparse.ArgumentParser(
  description='Generate figures from Monte Carlo data.',
  parents = [args_figs.parser])

args = parser.parse_args()

# directories, flags
datadir = './data/'
figdir = './figs/'
dataflag = '.dat'
figformat = '.pdf'

# figure options
max_RDF_radius = 10
dkte = args.ktemax / 1000
kte = arange(dkte,args.ktemax+dkte,dkte)

# load all data files
files = sort(glob.glob(datadir+'*'+dataflag))

# set figure parameters
xSize = 15
ySize = 10
fontSize = 12
lineWidth = 1
rcParams['font.size'] = fontSize
rcParams['font.size'] = fontSize
rcParams['axes.titlesize'] = fontSize
rcParams['legend.fontsize'] = fontSize
rcParams['legend.numpoints'] = 1
rcParams['xtick.labelsize'] = fontSize
rcParams['ytick.labelsize'] = fontSize
rcParams['lines.linewidth'] = lineWidth
rcParams['axes.color_cycle'] = ['k','b','g','r','c','m','y']

# figure template
def newFig():
    fig = figure(figsize=(xSize/2.54,ySize/2.54))
    ax = fig.add_subplot(1,1,1)
    fmt = mp.ticker.ScalarFormatter(useMathText=True,useOffset=False)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0,3))
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax

# interntal energy
def U(kte,energy,PDF):
    output = zeros(len(kte))
    for i in range(len(kte)):
        output[i] = sum(PDF*energy*exp(-energy/kte[i])) \
          / sum(exp(-energy/kte[i]))
    return output

# probability density and heat capacity plots
if (args.ww,args.ff,args.i) == (0,0,0) or args.all:
    # probability density over energy
    for file in [ file for file in files if '-E' in file ]:
        data = loadtxt(file)
        energy = -data[:,0][::-1]
        energy -= min(energy)
        PDF = data[:,1]/sum(data[:,1])
        fig, ax = newFig()
        plot(energy,PDF,'.')
        xlim(energy[0],energy[-1])
        xlabel('Energy')
        ylabel('Probability')
        tight_layout(pad=0.1)
        savename = figdir \
          + os.path.basename(file).replace(dataflag,figformat)
        savefig(savename)
        close()

    # heat capacity
    wws = []
    for file in files:
        ww = [ ww for ww in file.split('-') \
               if 'ww' in ww ][0].replace('ww','')
        if ww not in wws: wws.append(ww)
    for ww in wws:
        fig, ax = newFig()
        for file in [ file for file in files
                      if '-ww'+ww in file and '-E' in file ]:
            ids = os.path.basename(file).split('-')
            N = [ n for n in ids if 'N' in n][0].replace('N','')
            ff = [ ff for ff in ids if 'ff' in ff ][0].replace('ff','')
            data = loadtxt(file)
            energy = -data[:,0][::-1]
            energy -= min(energy)
            PDF = data[:,1]/sum(data[:,1])
            cv = (U(kte+dkte,energy,PDF) - U(kte-dkte,energy,PDF)) \
              / (2 * dkte) / int(N)
            plot(kte,cv,label='$\eta='+ff+'$')
            savename = '-'.join([ id for id in file.split('-')
                                  if 'ff' not in id ][:-1])
        savename = figdir+os.path.basename(savename)+'-cv'+figformat
        xlabel('$kT/\epsilon$')
        ylabel('$c_V/k$')
        xlim(0,args.ktemax+dkte)
        legend(loc=0)
        tight_layout(pad=0.1)
        savefig(savename)
        close()

# RDF plots
if (args.ww != 0 and args.ff != 0 and args.i != 0):
    for file in [ file for file in files if '-g' in file ]:
        data = loadtxt(file)
        energies = data[:,0]
        gs = data[:,1:]
        with open(file,'r') as stream:
            first_line = stream.readline().split(' ')
        for i in range(len(first_line)):
            if 'de_g' in first_line[i]:
                de_g = float(first_line[i+1])
                break
        radius = ((array(range(0,len(gs[0,:])))+0.5)*de_g)/2
        for i in range(len(energies)):
            if energies[i] >= args.i:
                plot(radius,gs[i,:],'.')
                axvline(1,color='k',linestyle=':')
                axvline(args.ww,color='k',linestyle=':')
                xlim(0,max_RDF_radius/2)
                xlabel('$r/\\sigma$')
                ylabel('$g(r)$')
                tight_layout(pad=0.1)
                savefig(figdir \
                        + name.replace(dataflag,'-i'+str(args.i)) \
                        + figformat)
                if args.show: show()
                close()
                exit(0)
