import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import rules

MCEXTRA = ['utilities.o', 'vector3d.o']

for o in MCEXTRA:
    rules.cxx(o[:-1]+'cpp')

for kind in ['square-well', 'polyhedra']:
    rules.cxx('Monte-Carlo/'+kind+'-monte-carlo.cpp')
    rules.cxx('Monte-Carlo/'+kind+'.cpp')
    rules.link('../'+kind+'-monte-carlo', ['Monte-Carlo/'+kind+'-monte-carlo.o',
                                     'Monte-Carlo/'+kind+'.o']+MCEXTRA)

for mc in ['soft-', 'pair-', 'triplet-', 'radial-distribution-', '']:
    rules.cxx('Monte-Carlo/'+mc+'monte-carlo.cpp')
    rules.link('../'+mc+'monte-carlo', ['Monte-Carlo/'+mc+'monte-carlo.o']+MCEXTRA)
