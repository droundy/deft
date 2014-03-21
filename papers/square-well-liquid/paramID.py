#!/usr/bin/python3
import os

# useful directories
thisdir = os.path.dirname(os.path.realpath(__file__))
jobdir = thisdir+'/jobs'

# Assumes this script is placed in [deft]/papers/square-well-liquid/
projectdir = os.path.realpath(thisdir+'../../..')
simname = 'square-well-monte-carlo'


# parameter ID class
class paramID:
    def __init__(self,walls,ww,ff,N,weights,files=[]):
        if walls == 0: self.walls = 'periodic'
        elif walls == 1: self.walls = 'walls'
        elif walls == 2: self.walls = 'tube'
        elif walls == 3: self.walls = 'box'
        self.ww = ww
        self.ff = ff
        self.N = N
        self.weights = weights
        for file in [ file for file in files
                      if self.name('N') in file ]:
            if 'E' in file: self.Efile = file
            if 'g' in file: self.gfile = file
    def __repr__(self):
        return str((self.walls,self.ww,self.ff,self.N))
    def name(self,option=''):
        out = self.walls + '-ww' + '%03.1f'%float(self.ww)
        if option:
            out += '-ff' + '%04.2f'%float(self.ff)
        if option == 'N':
            out += '-N' + str(self.N)
        return out + ('-nw' if not self.weights else '')
