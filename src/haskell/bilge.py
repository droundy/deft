#!/usr/bin/python2

import glob, re, string

importre = re.compile('^import\s+(\S+)', re.MULTILINE)
mainre = re.compile('^main\s+=', re.MULTILINE)

hsfiles = glob.glob('*.hs')

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
