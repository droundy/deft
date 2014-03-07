#!/usr/bin/python3

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
    def name(self,option='',interactions=None):
        out = self.walls + '-ww' + '%03.1f'%float(self.ww)
        if option:
            out += '-ff' + '%04.2f'%float(self.ff)
        if option == 'N':
            out += '-N' + str(self.N)
            if interactions != None:
                out += '-i' + str(interactions)
        return out + ('-nw' if not self.weights else '')
