#!/usr/bin/python3

import config

import os, re, sys

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
    def compile(self, cpp):
        obj = cpp[:-3]+'o'
        print('\n|', config.cxx, config.cxxflags, '-c', '-o', obj, cpp, file=self._f)
        print('<', cpp, file=self._f)
        print('>', obj, file=self._f)
    def link(self, maincpp, exe):
        obj = {maincpp[:-3]+'o'}
        with open(maincpp) as f:
            maincontents = f.read()
            for i in includes_re.findall(maincontents):
                if i == 'Monte-Carlo/monte-carlo':
                    i = 'utilities'
                print('# include ', i, file=self._f)
                if os.path.exists('src/'+i+'.cpp'):
                    print('# should link with', 'src/'+i+'.cpp', file=self._f)
                    obj |= {'src/'+i+'.o'}
                else:
                    print('# should not link with', 'src/'+i+'.cpp', file=self._f)
            print('# instructions:', instructions_re.findall(maincontents), file=self._f)
            for i in instructions_re.findall(maincontents):
                print('# instructions:', repr(i), file=self._f)
                try:
                    exec(i)
                except:
                    print('# instructions failed!', repr(i), sys.exc_info(), file=self._f)
        print('\n|', config.cxx, config.linkflags, '-o', exe, ' '.join(obj), file=self._f)
        print('>', exe, file=self._f)
        for o in obj:
            print('<', o, file=self._f)
