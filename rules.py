import os.path as path
import string

CXXFLAGS='-g -O3 -std=c++11 -I . -I ../include'

generic_sources = ['src/lattice.o', 'src/utilities.o',
                   'src/Faddeeva.o', 'src/GridDescription.o',
                   'src/Grid.o', 'src/ReciprocalGrid.o',
                   'src/IdealGas.o', 'src/ChemicalPotential.o',
                   'src/HardSpheres.o', 'src/ExternalPotential.o',
                   'src/Functional.o', 'src/ContactDensity.o',
                   'src/Gaussian.o', 'src/Pow.o',
                   'src/WaterSaftFast.o',
                   'src/WaterSaft_by_handFast.o',
                   'src/EffectivePotentialToDensity.o',
                   'src/equation-of-state.o', 'src/water-constants.o',
                   'src/compute-surface-tension.o',
                   'src/Minimizer.o', 'src/Downhill.o',
                   'src/Precision.o', 'src/ConjugateGradient.o',
                   'src/QuadraticLineMinimizer.o',
                   'src/SteepestDescent.o']

def cxx(cpp):
    #print '| pwd'
    obj = cpp[:-3]+'o'
    if cpp[-3:] == '.cc':
        obj = cpp[:-2]+'o'
    print '# cpp', cpp, obj, cpp[-3:]
    print '| g++ ' + CXXFLAGS + ' -c %s -o %s' % (cpp, obj)
    print '>', obj
    print

def link(binary, objects):
    print '| g++ -lfftw3 -lpopt -o %s %s' % (binary, string.join(objects))
    print '> %s' % binary
    for o in objects:
        print '<', o
    print

def linkgeneric(binary, objects):
    link(binary, objects + generic_sources)

def latex(tex):
    d = path.dirname(tex)
    t = path.basename(tex)
    b = t[:-4]
    print ('| cd %s && pdflatex %s && bibtex %s && pdflatex %s && pdflatex %s'
            % (d, t, b, t, t))
    print "> %s.pdf" % b
    print

def ghc_hs(module):
    print "| ghc -O2 -c -o %s.o %s.hs" % (module, module)
    print "> %s.o" % module
    print "> %s.hi" % module
    print

def ghclink(binary, source, objects):
    print "| ghc -O2 --make -o %s %s %s" % (binary, source, string.join(objects))
    print "> %s.o" % module
    for ob in objects:
        print '< %s', o
    print

def pyc(py):
    print "| python -m compileall %s" % py
    print "> %sc" % py
    print
