#!/usr/bin/python

import glob, re, string

importre = re.compile('^import\s+(\S+)', re.MULTILINE)
mainre = re.compile('^main\s+=', re.MULTILINE)

hsfiles = glob.glob('*.hs')

imports = {}
mainfiles = []

allobjects = [hsf[:-2]+'o' for hsf in hsfiles]
allhi = [hsf[:-2]+'hi' for hsf in hsfiles]

for hsf in hsfiles:
    f = open(hsf, 'r')
    hs = f.read()
    f.close()
    if mainre.search(hs):
        mainfiles.append(hsf)
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
        if isokay:
            print ': %s | %s |> !ghc |>' % (hsf, string.join([i + '.hi' for i in imports[hsf]]))
            hsfiles.remove(hsf)

allobjects = string.join(allobjects + allhi)
for main in mainfiles:
    print ': %s | %s |> !ghclink |>' % (main, allobjects)
