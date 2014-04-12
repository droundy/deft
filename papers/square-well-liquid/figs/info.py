#!/usr/bin/python2
import os, glob
from pylab import sort
import matplotlib as mp

####################
# file info
####################

# directories, flags, extensions
figdir = os.path.dirname(os.path.realpath(__file__))
swdir = os.path.dirname(figdir)
projectdir = os.path.realpath(swdir+'../../..')
jobdir = swdir+'/jobs'
datadir = swdir+'/data/'
simname = 'square-well-monte-carlo'
data_ext = '.dat'
fig_ext = '.pdf'

for dir in [jobdir,datadir,figdir]:
    if not os.path.isdir(dir):
        os.mkdir(dir)

# load all data files
files = sort(glob.glob(datadir+'*'+data_ext))
E_files = [ file for file in files if '-E' in file ]
g_files = [ file for file in files if '-g' in file ]

wall_options = [(0,'periodic'),(1,'walls'),(2,'tube'),(3,'box')]
def wall_swap(input): # swap wall number and tag
    for option in wall_options:
        if input == option[0]: return option[1] # number to tag
        if option[1] in input: return option[0] # tag to number
    return 'error'

# parameter ID class
class paramID:
    def __init__(self,walls,ww,ff,N,nw,files=[]):
        self.walls = wall_swap(walls)
        self.ww = ww
        self.ff = ff
        self.N = N
        self.nw = nw # nw = "no weights"
        for file in [ file for file in files
                      if self.name('ff-N') in file ]:
            if 'E' in file: self.Efile = file
            if 'g' in file: self.gfile = file
    def __repr__(self):
        return str((self.walls,self.ww,self.ff,self.N))
    def name(self,option=''):
        out = self.walls + '-ww' + '%04.2f'%float(self.ww)
        if 'ff' in option:
            out += '-ff' + '%04.2f'%float(self.ff)
        if 'N' in option:
            out += '-N' + str(self.N)
        return out + ('-nw' if self.nw else '')

# store possible parameters and create list of sets to plot
wws = []
ffs = []
Ns = []
param_list = []
for file in E_files:
    ids = file.split('-')
    for id in ids:
        if 'ww' in id: ww = float(id.replace('ww',''))
        if 'ff' in id: ff = float(id.replace('ff',''))
        if 'N' in id: N = int(id.replace('N',''))
    walls = wall_swap(os.path.basename(file))
    nw = ('nw' in os.path.basename(file))
    if ww not in wws:
        wws.append(ww)
    if ff not in ffs:
        ffs.append(ff)
    if N not in Ns:
        Ns.append(N)
    p = paramID(walls,ww,ff,N,nw,files)
    if p not in param_list: param_list.append(p)

####################
# figure options
####################

mp.rcParams['legend.numpoints'] = 3
mp.rcParams['lines.linewidth'] = 1
mp.rcParams['lines.markersize'] = 2
mp.rcParams['axes.color_cycle'] = ['k','b','g','r','c','m','y']

default_font_size = 12 # pt
def use_font(font_size):
    mp.rcParams['font.size'] = font_size
    mp.rcParams['axes.titlesize'] = font_size
    mp.rcParams['legend.fontsize'] = font_size
    mp.rcParams['xtick.labelsize'] = font_size
    mp.rcParams['ytick.labelsize'] = font_size
use_font(default_font_size)

# figure template
default_size = (8, 5) # (x,y) in cm
default_power_limit = (-2, 3)
def new_fig(size=default_size,power_limit=default_power_limit,
           use_power_limits=True,use_offset=False):
    fig = mp.pyplot.figure(figsize=size)
    ax = fig.add_subplot(1,1,1)
    fmt = mp.ticker.ScalarFormatter(useMathText=True,useOffset=use_offset)
    if use_power_limits:
        fmt.set_powerlimits(power_limit)
    fmt.set_scientific(True)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax
