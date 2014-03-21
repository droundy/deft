#!/usr/bin/python3
import matplotlib as mp
mp.use('Agg')
from pylab import *
import os, sys, glob, socket, argparse

import arguments
from paramID import *

# parse arguments
parser = argparse.ArgumentParser(
  description='Generate figures from Monte Carlo data.',
  parents = [arguments.parser])

parser.add_argument(
		'-ktemin', metavar='FLOAT', type=float, default=0,
    help='Maximum kt/e in heat capacity plot')

parser.add_argument(
		'-ktemax', metavar='FLOAT', type=float, default=30,
    help='Maximum kt/e in heat capacity plot')

parser.add_argument(
		'-kteplot', metavar='INT', type=int, nargs='+',
    default=[3,10,30],
    help='Value(s) of kt/e to plot for radial distribution function')

parser.add_argument(
		'--print', action='store_true', help='Print more status text')

parser.add_argument(
		'--dos', action='store_true',
    help='Generate density of states figures')

parser.add_argument(
		'--hc', action='store_true',
    help='Generate heat capacity figures')

parser.add_argument(
		'--rd', action='store_true',
    help='Generate radial distribution figures')

args = parser.parse_args()

# directories, flags
datadir = './data/'
figdir = './figs/'
dataflag = '.dat'
figformat = '.pdf'

# figure options
max_RDF_radius = 10
dkte = (args.ktemax-args.ktemin) / 1000
kte = arange(args.ktemin+dkte,args.ktemax,dkte)

# set figure parameters
xSize = 15 # cm
ySize = 10 # cm
fontSize = 12 # pt
rcParams['font.size'] = fontSize
rcParams['font.size'] = fontSize
rcParams['axes.titlesize'] = fontSize
rcParams['legend.fontsize'] = fontSize
rcParams['legend.numpoints'] = 3
rcParams['xtick.labelsize'] = fontSize
rcParams['ytick.labelsize'] = fontSize
rcParams['lines.linewidth'] = 1
rcParams['lines.markersize'] = 2
rcParams['axes.color_cycle'] = ['k','b','g','r','c','m','y']

# load all data files
files = sort(glob.glob(datadir+'*'+dataflag))

# deside which figures to generate
if not args.dos and not args.hc and not args.rd:
    args.dos = True
    args.hc = True
    args.rd = True

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
    p = paramID(args.walls,ww,ff,N,args.nw,files)
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
    fmt.set_powerlimits((-2,3))
    fmt.set_scientific(True)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax

# interntal energy relative to well depth
def Ue(kte,counts,DS):
    output = zeros(len(kte))
    for i in range(len(output)):
        output[i] = sum(counts*DS*exp(-counts/kte[i])) \
          / sum(DS*exp(-counts/kte[i]))
    return output

# probability density over energy
if args.dos:
    print("Generating density of states figures")
    for ww in args.ww:
        if args.print: print('  well width:',ww)
        fig, ax = newFig()
        for p in [ p for p in paramList
                   if p.ww == ww and p.ff in args.ff ]:
                data = loadtxt(p.Efile,ndmin=2)
                energy = -data[:,0][::-1]/p.N
                DS = data[:,1]
                DS /= sum(DS)
                p.max_index = argmax(DS)
                semilogy(energy,DS,'.',label=figLabel(p,'N'))
        xlabel('$E/N\epsilon$')
        ylabel('$D$')
        legend(loc=0)
        tight_layout(pad=0.1)
        savefig(figdir+p.name()+'-dos'+figformat)
        close()

# heat capacity
if args.hc:
    print("Generating heat capacity figures")
    for ww in args.ww:
        if args.print: print('  well width:',ww)
        fig, ax = newFig()
        for p in [ p for p in paramList
                   if p.ww == ww and p.ff in args.ff ]:
                data = loadtxt(p.Efile,ndmin=2)
                counts = data[:,0][::-1]
                counts -= min(counts)
                DS = data[:,1]
                cv = (Ue(kte+dkte/2,counts,DS)
                      - Ue(kte-dkte/2,counts,DS)) / dkte / p.N
                plot(kte,cv,label=figLabel(p,'N'))
        xlabel('$kT/\epsilon$')
        ylabel('$C_V/Nk$')
        xlim(args.ktemin,args.ktemax)
        legend(loc=0)
        tight_layout(pad=0.1)
        savefig(figdir+p.name()+'-hc'+figformat)
        close()

# radial distribution function at fixed temperature
def g(gs,kte,counts,DS):
    output = zeros(len(gs[0,:]))
    Z = sum(DS*exp(-counts/kte))
    for i in range(len(output)):
        output[i] += sum(gs[:,i]*DS*exp(-counts/kte)) / Z
    return output

# radial distribution function
if args.rd:
    print("Generating radial distribution figures")
    for ww in args.ww:
        for p in [ p for p in paramList
                   if p.ww == ww and p.ff in args.ff ]:
            if args.print: print('  trial:',p.name('N'))
            g_data = loadtxt(p.gfile,ndmin=2)
            g_counts = g_data[:,0][::-1]
            gs = g_data[:,1:]

            E_data = loadtxt(p.Efile,ndmin=2)
            counts = E_data[:,0][::-1]
            use_counts = [ i for i in range(len(counts))
                           if counts[i] in g_counts ]
            counts = counts[use_counts]
            counts -= min(counts)
            DS = E_data[use_counts,1]
            DS /= sum(DS)

            with open(p.gfile,'r') as stream:
                first_line = stream.readline().split(' ')
            for i in range(len(first_line)):
                if 'de_g' in first_line[i]:
                    de_g = float(first_line[i+1])
                    break
            radius = (array(range(0,len(gs[0,:])))+0.5) * de_g/2

            fig, ax = newFig()
            for kte in args.kteplot:
                plot(radius,g(gs,kte,counts,DS),'.',
                     label='$kT/\epsilon='+str(kte)+'$')
            axvline(1,color='k',linestyle=':')
            axvline(p.ww,color='k',linestyle=':')
            xlim(0,max_RDF_radius/2)
            xlabel('$r/\\sigma$')
            ylabel('$g(r)$')
            legend(loc=0)
            tight_layout(pad=0.1)
            savefig(figdir+p.name('N')+'-rd'+figformat)
            close()
