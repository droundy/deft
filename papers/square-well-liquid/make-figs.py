#!/usr/bin/python3
from __future__ import division
import matplotlib as mp
#mp.use('Agg')
from pylab import *
import os, sys, glob, socket, argparse

import figArgs

# parse arguments
parser = argparse.ArgumentParser(
  description='Generate figures from Monte Carlo data.',
  parents = [figArgs.parser])

args = parser.parse_args()

# directories, flags
datadir = './data/'
figdir = './figs/'
dataflag = '.dat'
figformat = '.pdf'

# figure options
max_RDF_radius = 10

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
    fmt = mp.ticker.ScalarFormatter(useMathText=True)
    fmt.set_scientific(True)
    fmt.set_powerlimits((0,3))
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax

# single simulation figures
for file in files:
    name = os.path.basename(file)

    if (args.ff,args.ww,args.i) == (0,0,0) and '-E' in file:
        # energy probability density plots
        data = loadtxt(file)
        energy = -data[:,0][::-1]
        PDF = data[:,1]/sum(data[:,1])
        fig, ax = newFig()
        plot(energy,PDF,'.')
        xlim(energy[0],energy[-1])
        xlabel('Energy')
        ylabel('Probability')
        tight_layout(pad=0.1)
        savefig(figdir+name.replace(dataflag,figformat))
        close()

    elif (args.ff,args.ww,args.i) != (0,0,0) \
      and '-g' in file \
        and 'ff'+('%4.2f'%args.ff) in file \
        and 'ww'+('%i'%args.ww) in file:
        # RDF plots
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
                savefig(figdir
                        + name.replace(dataflag,'-i'+str(args.i))
                        + figformat)
                if args.show:
                    show()
                close()
                exit(0)
