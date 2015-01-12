import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

import glob
import rules

MCEXTRA = ['utilities.o', 'vector3d.o']

os.chdir('src')
for cpp in glob.glob('*.cpp') + glob.glob('*.cc'):
    rules.cxx(cpp)

for kind in ['square-well', 'polyhedra']:
    rules.cxx('Monte-Carlo/'+kind+'-monte-carlo.cpp')
    rules.cxx('Monte-Carlo/'+kind+'.cpp')
    rules.link('../'+kind+'-monte-carlo', ['Monte-Carlo/'+kind+'-monte-carlo.o',
                                     'Monte-Carlo/'+kind+'.o']+MCEXTRA)

for mc in ['soft-', 'pair-', 'triplet-', 'radial-distribution-', '']:
    rules.cxx('Monte-Carlo/'+mc+'monte-carlo.cpp')
    rules.link('../'+mc+'monte-carlo', ['Monte-Carlo/'+mc+'monte-carlo.o']+MCEXTRA)

for x in ['HardSpheresNoTensor2Fast',
          'TensorWhiteBearFast',
          'WhiteBearMarkIIFast',
          'gSigmaA_by_handFast',
          'gSigmaA_automagicFast',
          'TensorDensityXXFast',
          'n2DensityFast',
          'VectorDensityXFast',
          'YuWuCorrelationFast',
          'SaftFluid2Fast',
          'EntropySaftFluid2Fast',
          'CorrelationGrossCorrectFast',
          'gSigmaSm2Fast',
          'gSigmaAm2Fast',
          'gSigmaSFast',
          'gSigmaAFast',
          'HughesHBFast',
          'SoftFluidFast',
          'HardFluidFast',
          'HardRosenfeldFluidFast',
          'WaterXFast',
          'HughesXFast']:
    print '| cd .. && src/haskell/functionals.exe src/%s' % (x + '.cpp')
    print '< haskell/functionals.exe'
    print '> %s' % (x + '.cpp')

