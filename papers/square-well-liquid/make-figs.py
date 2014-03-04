#!/usr/bin/python3
from __future__ import division
import matplotlib as mp
#mp.use('Agg')
from pylab import *
import os, glob, socket

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
    data = loadtxt(file,ndmin=2)
    fig, ax = newFig()

    # energy histograms
    if '-E-' in file:
        plot(-data[:,0],data[:,1]/sum(data[:,1]),'.')
        xlim(-data[-1,0],-data[1,0])
        xlabel('Energy')
        ylabel('Probability')
        tight_layout(pad=0.1)
        savefig(figdir+name.replace(dataflag,figformat))
        close()

    elif '-g-' in file:
        radii = len(data[0,:])-1
        with open(file,'r') as f:
            first_line = f.readline().split(' ')
        for i in range(len(first_line)):
            if 'de_g' in first_line[i]:
                de_g = float(first_line[i+1])
                break
        radius = (array(range(0,radii))+0.5)*de_g
        for i in range(0, len(data[:,0]), 10):
            plot(radius,data[i,1:])
            xlim(0,max_RDF_radius)
            xlabel('$r/R$')
            ylabel('$g(r)$')
            tight_layout(pad=0.1)
            savefig(figdir \
                    + name.replace(dataflag,
                                   '-i-'+str(int(data[i,0]))+figformat))
            close()
