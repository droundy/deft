#!/usr/bin/python3

import config

import os, re, sys
import numpy
import numpy as np

includes_re = re.compile(r'#include\s+"(\S+)\.h"')

instructions_re = re.compile(r'/\*--(.+)--\*/', re.MULTILINE | re.DOTALL)

while not os.path.exists('.git'):
    os.chdir('..') # go back into root directory of source code

class facfile:
    def __init__(self, f):
        self._facname = f
        self._f = open(f, 'w')
        self._dirname = os.path.dirname(f)
    def rule(self, cmd, inputs, outputs):
        print('\n?', cmd, file=self._f)
        for i in inputs:
            print('<', i, file=self._f)
        for o in outputs:
            print('>', o, file=self._f)
    def default(self, cmd, inputs, outputs):
        print('\n|', cmd, file=self._f)
        for i in inputs:
            print('<', i, file=self._f)
        for o in outputs:
            print('>', o, file=self._f)
    def compile(self, cpp):
        obj = cpp[:-3]+'o'
        print('\n|', config.cxx, config.cxxflags, '-c', '-o', obj, cpp, file=self._f)
        print('<', cpp, file=self._f)
        print('>', obj, file=self._f)
    def link(self, maincpp, exe, objects=set([])):
        obj = {maincpp[:-3]+'o'}
        headers = set([])
        with open(maincpp) as f:
            maincontents = f.read()
            for i in includes_re.findall(maincontents):
                if i == 'Monte-Carlo/monte-carlo':
                    i = 'utilities'
                if i[-4:] == 'Fast' or os.path.exists('src/'+i+'.cpp'):
                    obj |= {'src/'+i+'.o'}
                    headers |= {'src/'+i+'.h'}
            for i in instructions_re.findall(maincontents):
                try:
                    exec(i)
                except:
                    print('# instructions failed!', repr(i), file=self._f)
                    print('# instructions failed!', sys.exc_info()[1], file=self._f)
        print('\n|', config.cxx, config.linkflags, '-o', exe, ' '.join(obj | objects), file=self._f)
        print('>', exe, file=self._f)
        for o in obj | objects | headers:
            print('<', o, file=self._f)
