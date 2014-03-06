#!/usr/bin/python3
import matplotlib as mp
#mp.use('Agg')
from pylab import *
import os, sys, glob, socket, argparse

import arguments, paramID

# parse arguments
parser = argparse.ArgumentParser(
  description='Generate figures from Monte Carlo data.',
  parents = [arguments.parser])

parser.add_argument(
		'-i', metavar='INT', type=int, nargs='+', default=[],
    help='Number(s) of interactions')

parser.add_argument(
		'-ktemax', metavar='FLOAT', type=float, default=150,
    help='Maximum kt in heat capacity plot')

parser.add_argument(
		'--show', action='store_true',
    help='Show generated RDF figure')

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

# load all data files
files = sort(glob.glob(datadir+'*'+dataflag))
rcParams['axes.color_cycle'] = ['k','b','g','r','c','m','y']

# store possible parameters and create list of sets to plot
wws = []
ffs = []
Ns = []
paramList = []
for file in [ file for file in files if '-E' in file ]:
    ids = file.split('-')
    for id in ids:
        if 'ww' in id: ww = float(id.replace('ww',''))
        if 'ff' in id: ff = float(id.replace('ff',''))
        if 'N' in id: N = int(id.replace('N',''))
    if ww not in wws:
        wws.append(ww)
    if ff not in ffs:
        ffs.append(ff)
    if N not in Ns:
        Ns.append(N)
    p = paramID.paramID(args.walls,ww,ff,N,args.weights,files)
    if p not in paramList: paramList.append(p)

# if an parameter isn't given, assume all possible values are wanted
if args.ww == []: args.ww = wws
if args.ff == []: args.ff = ffs
if args.N == []: args.N = Ns

# figure label format
def figLabel(p,option=''):
    return '$\eta=' + str(p.ff) + '$' \
      + (', $N=' + str(p.N) + '$' if option == 'N' else '')

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

# interntal energy relative to well depth
def Ue(kte,counts,PDF):
    output = zeros(len(kte))
    for i in range(len(kte)):
        if kte[i] != 0:
            output[i] = sum(counts*PDF*exp(-counts/kte[i])) \
              / sum(PDF*exp(-counts/kte[i]))
    return output

if args.i == []:
    # probability density over energy
    print("Generating probability distribution plots")
    for ww in args.ww:
        fig, ax = newFig()
        for p in [ p for p in paramList
                   if p.ww == ww
                   and p.ff in args.ff ]:
            data = loadtxt(p.Efile,ndmin=2)
            energy = -data[:,0][::-1]
            energy -= min(energy)
            PDF = data[:,1]/sum(data[:,1])
            plot(energy,PDF,'.',label=figLabel(p,'N'))
        xlabel('Energy ($\epsilon/N$)')
        ylabel('Probability')
        legend(loc=0)
        tight_layout(pad=0.1)
        savefig(figdir+p.name()+'-PDF'+figformat)
        close()

    # heat capacity
    print("Generating heat capacity plots")
    for ww in args.ww:
        fig, ax = newFig()
        for p in [ p for p in paramList
                   if p.ww == ww
                   and p.ff in args.ff ]:
            data = loadtxt(p.Efile,ndmin=2)
            counts = data[:,0][::-1]
            counts -= min(counts)
            PDF = data[:,1]/sum(data[:,1])
            cv = (Ue(kte+dkte,counts,PDF) - Ue(kte-dkte,counts,PDF)) \
              / (2 * dkte) / int(N)
            plot(kte,cv,label=figLabel(p,'N'))
        xlabel('$kT/\epsilon$')
        ylabel('$c_v/k$')
        xlim(0,args.ktemax+dkte)
        legend(loc=0)
        tight_layout(pad=0.1)
        savefig(figdir+p.name()+'-cv'+figformat)
        close()

# radial distribution function
for energy in args.i:
    print("Generating radial distribution plots for",
          energy,"interactions")
    for ww in args.ww:
        for ff in args.ff:
            for N in args.N:
                for p in [ p for p in paramList
                           if p.ww == ww and p.ff == ff and p.N == N ]:
                    data = loadtxt(p.gfile)
                    energies = data[:,0]
                    gs = data[:,1:]
                    with open(p.gfile,'r') as stream:
                        first_line = stream.readline().split(' ')
                    for i in range(len(first_line)):
                        if 'de_g' in first_line[i]:
                            de_g = float(first_line[i+1])
                            break
                    radius = (array(range(0,len(gs[0,:])))+0.5) \
                              * de_g/2
                    for i in range(len(energies)):
                        if energies[i] >= energy:
                            plot(radius,gs[i,:],'.')
                            break
                    axvline(1,color='k',linestyle=':')
                    axvline(ww,color='k',linestyle=':')
                    xlim(0,max_RDF_radius/2)
                    xlabel('$r/\\sigma$')
                    ylabel('$g(r)$')
                    tight_layout(pad=0.1)
                    savefig(figdir+p.name('N',energy)+figformat)
                    close()
