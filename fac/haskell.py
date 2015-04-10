#!/usr/bin/python3

import re, string, os

import facfile

haskell = facfile.facfile('src/haskell/.haskell.fac')

os.chdir('src/haskell')

importre = re.compile('^import\s+(\S+)', re.MULTILINE)
mainre = re.compile('^main\s+=', re.MULTILINE)

hsfiles = """
  C CodeGen Expression ExternalLennardJones ExternalPotentialTest
  FMT functionals HughesSaft IdealGas LatexDouble latex-functionals
  Latex LogN0 NewCode Optimize QuadraticGaussian Quadratic QuadraticN0
  Rosenfeld SFMT Statement test WaterSaft WhiteBear SW_liquid
""".split()

all_objects = [x+'.o' for x in hsfiles]

generated_names = """
  ExternalPotentialTest HomogeneousSFMTFluid HomogeneousWaterSaftByHand
  HomogeneousWaterSaft HomogeneousWhiteBear
  LogN0 Phi1 Phi2 Phi3 QuadraticGaussian Quadratic QuadraticN0
  SFMTFluid SFMTFluidVeff
  SPhi1 SPhi2 SPhi3
  SW_liquid
  WaterSaftByHand WaterSaft WhiteBearFluid WhiteBear
""".split()

old_generated = """
  HardSpheresNoTensor2Fast TensorWhiteBearFast WhiteBearMarkIIFast
  gSigmaA_by_handFast gSigmaA_automagicFast TensorDensityXXFast
  n2DensityFast VectorDensityXFast YuWuCorrelationFast SaftFluid2Fast
  EntropySaftFluid2Fast CorrelationGrossCorrectFast gSigmaSm2Fast
  gSigmaAm2Fast gSigmaSFast gSigmaAFast HughesHBFast SoftFluidFast
  HardFluidFast HardRosenfeldFluidFast WaterXFast HughesXFast
""".split()


imports = {}
mainfiles = []

for hsf in hsfiles:
    f = open(hsf+'.hs', 'r')
    hs = f.read()
    f.close()
    if mainre.search(hs):
        mainfiles.append(hsf)
    imports[hsf] = importre.findall(hs)
    imports[hsf] = [i for i in imports[hsf] if i in hsfiles]

done = {}
while len(hsfiles) > 0:
    for hsf in hsfiles:
        isokay = True
        for i in imports[hsf]:
            if i+'.hs' in hsfiles:
                isokay = False
                break
        if isokay:
            haskell.rule('ghc -O2 -c %s.hs' % hsf,
                         [x+'.hi' for x in imports[hsf]],
                         [hsf+'.o', hsf+'.hi'])
            hsfiles.remove(hsf)

def all_dependencies(x):
    handled = {x}
    deps = {x} | set(imports[x])
    while len(deps) > len(handled):
        for d in deps - handled:
            deps |= set(imports[d])
            handled |= {d}
    return deps

for main in mainfiles:
    objects = [o+'.o' for o in all_dependencies(main)]
    haskell.rule('ghc -O2 -package containers -package filepath -package directory -o %s.exe %s'
                 % (main, ' '.join(objects)), objects, [main+'.exe'])

haskell.rule('python2 create_generators.py',
             [], ['generate_%s.hs' % x for x in generated_names])

os.chdir('../..')
cxx = facfile.facfile('.generated-code.fac')

for name in generated_names:
    # command to compile haskell code
    haskell.rule('ghc -O2 --make -o generate_%s.exe generate_%s.hs' % (name, name),
                 ['generate_%s.hs' % name] + all_objects,
                 ['generate_%s.o' % name, 'generate_%s.hi' % name,
                  'generate_%s.exe' % name])

    # the following "if" avoid automatically regenerating cpp code
    # that takes a long long time to generate.
    if name not in ['SW_liquid']:
        # command to generate C++ code
        haskell.rule('cd ../.. && src/haskell/generate_%s.exe' % name,
                     ['generate_%s.exe' % name],
                     ['../new/%sFast.cpp' % name, '../new/%sFast.h' % name])

    # command to compile C++ code
    cxx.compile('src/new/%sFast.cpp' % name)

for name in old_generated:
    # command to generate C++ code
    haskell.rule('cd ../.. && src/haskell/functionals.exe src/%s.cpp' % name,
                 ['functionals.exe'],
                 ['../%s.cpp' % name])
    # command to compile C++ code
    cxx.compile('src/%s.cpp' % name)

