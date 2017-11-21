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
    def rule(self, cmd, inputs, outputs, Caches=['/etc/ld.so.cache'], caches=[]):
        print('\n?', cmd, file=self._f)
        for i in inputs:
            print('<', i, file=self._f)
        for o in outputs:
            print('>', o, file=self._f)
        for C in Caches:
            print('C', C, file=self._f)
        for c in caches:
            print('c', c, file=self._f)
    def default(self, cmd, inputs, outputs, Caches=['/etc/ld.so.cache'], caches=[]):
        print('\n|', cmd, file=self._f)
        for i in inputs:
            print('<', i, file=self._f)
        for o in outputs:
            print('>', o, file=self._f)
        for C in Caches:
            print('C', C, file=self._f)
        for c in caches:
            print('c', c, file=self._f)
    def compile(self, cpp):
        obj = cpp[:-3]+'o'
        print('\n|', config.cxx, config.cxxflags, '-c', '-o', obj, cpp, file=self._f)
        print('<', cpp, file=self._f)
        print('>', obj, file=self._f)
    def compile_big_file(self, cpp):
        obj = cpp[:-3]+'o'
        slow_flags = ['-flto', '-O3', '-g']
        add_flags = ['-O1']
        print('\n|', config.cxx,
              ' '.join(add_flags+
                       [f for f in config.cxxflags.split()
                        if f not in slow_flags]),
              '-c', '-o', obj, cpp, file=self._f)
        print('<', cpp, file=self._f)
        print('>', obj, file=self._f)
    def link(self, maincpp, exe, objects=set([]), extra_requirements=set([])):
        obj = {maincpp[:-3]+'o'}
        headers = set([])
        with open(maincpp) as f:
            maincontents = f.read()
            includes = includes_re.findall(maincontents)
            if "OptimizedFunctionals" in includes:
                obj |= set(""" src/SoftFluidFast.o src/SaftFluid2Fast.o
                               src/EntropySaftFluid2Fast.o src/WaterSaftFast.o src/WaterXFast.o
                               src/HughesXFast.o src/HughesHBFast.o src/TensorWhiteBearFast.o
                               src/WhiteBearMarkIIFast.o src/HardSpheresNoTensor2Fast.o
                               src/YuWuCorrelationFast.o
                            """.split())
            if "ContactDensity" in includes:
                obj |= set(""" src/gSigmaAFast.o src/gSigmaSFast.o
                               src/gSigmaAm2Fast.o src/gSigmaSm2Fast.o
                               src/YuWuCorrelationFast.o
                            """.split())
            if "version-identifier" in includes:
                headers |= {'src/version-identifier.h'}
            for i in includes:
                if i == 'Monte-Carlo/monte-carlo':
                    i = 'utilities'
                if i[-4:] == 'Fast' or os.path.exists('src/'+i+'.cpp'):
                    obj |= {'src/'+i+'.o'}
                    headers |= {'src/'+i+'.h'}
            for i in instructions_re.findall(maincontents):
                try:
                    # Ideally we would like to add /etc/ld.so.cache as
                    # cache when running code.  However, this would
                    # require editing all our .cpp files.  :(
                    exec(i)
                except:
                    print('# instructions failed for', maincpp, repr(i), file=self._f)
                    print('# instructions failed!', sys.exc_info()[1], file=self._f)
        print('\n|', config.cxx, '-o', exe,
              ' '.join(sorted(obj | objects)), config.linkflags, file=self._f)
        print('>', exe, file=self._f)
        for o in obj | objects | headers:
            print('<', o, file=self._f)
        mainobj = maincpp[:-3]+'o'
        print('\n|', config.cxx, config.cxxflags, '-c', '-o', mainobj, maincpp, file=self._f)
        print('<', maincpp, file=self._f)
        print('>', mainobj, file=self._f)
        for h in headers | extra_requirements:
            print('<', h, file=self._f)
