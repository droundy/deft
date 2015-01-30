#!/usr/bin/python2

import glob, re, string

importre = re.compile('^import\s+(\S+)', re.MULTILINE)
mainre = re.compile('^main\s+=', re.MULTILINE)

hsfiles = [hs for hs in glob.glob('*.hs') if hs[:9] != 'generate_']

imports = {}
mainfiles = []

allobjects = [hsf[:-2]+'o' for hsf in hsfiles]

for hsf in hsfiles:
    f = open(hsf, 'r')
    hs = f.read()
    f.close()
    if mainre.search(hs):
        mainfiles.append(hsf)
        allobjects.remove(hsf[:-2]+'o')
    imports[hsf] = importre.findall(hs)
    imports[hsf] = [i for i in imports[hsf] if i+'.hs' in hsfiles]

done = {}
while len(hsfiles) > 0:
    for hsf in hsfiles:
        isokay = True
        for i in imports[hsf]:
            if i+'.hs' in hsfiles:
                isokay = False
                break
        if hsf in mainfiles:
            hsfiles.remove(hsf)
        elif isokay:
            print '| ghc -O2 -c %s' % hsf
            print '> %s.o' % hsf[:-3]
            print '> %s.hi' % hsf[:-3]
            for i in imports[hsf]:
                print '< %s.hi' % i
            print
            hsfiles.remove(hsf)

for main in mainfiles:
    exe = main[:-3]+'.exe'
    print '| ghc -O2 --make -o %s %s' % (exe, main)
    print '>', exe
    for o in allobjects:
        print '<', o
    print

generated_names = ["ExternalPotentialTest",
                   "Quadratic",
                   "QuadraticN0",
                   "QuadraticGaussian",
                   "LogN0",
                   "Phi1",
                   "Phi2",
                   "Phi3",
                   "SPhi1",
                   "SPhi2",
                   "SPhi3",
                   "HomogeneousWhiteBear",
                   "WhiteBear",
                   "WhiteBearFluid",
                   "SFMTFluid",
                   "SFMTFluidVeff",
                   "HomogeneousSFMTFluid",
                   "WaterSaft",
                   "WaterSaftByHand",
                   "HomogeneousWaterSaft",
                   "HomogeneousWaterSaftByHand"]

print '| python2 create_generators.py'
print '< create_generators.py'
for n in generated_names:
    print '> generate_%s.hs' % n

for name in generated_names:
    # command to compile haskell code
    print '| ghc -O2 --make -o generate_%s.exe generate_%s.hs' % (name, name)
    print '> generate_%s.exe' % name
    print '< generate_%s.hs' % name
    for o in allobjects:
        print '<', o
    print
    # command to generate C++ code
    print '| cd ../.. && src/haskell/generate_%s.exe' % name
    print '> ../new/%sFast.cpp' % name
    print '> ../new/%sFast.h' % name
    print '< generate_%s.exe' % name
    print
