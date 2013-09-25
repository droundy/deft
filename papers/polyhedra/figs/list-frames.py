#!/usr/bin/python
import glob

names = glob.glob("figs/mc/vertices/*--1.dat")
bases = [name[:-6] for name in names]

info = []
for base in bases:
  frames = 0
  for i in glob.glob(base+"*"):
    if len(i.split('-')) == 6: #number is positive
      frames += 1
  name = base.split('-')
  celltype = name[0].split('/')[-1]
  ff = name[1]
  polyhedron = name[3]
  N = name[4]
  info.append(dict(celltype=celltype, name=polyhedron, ff=ff, N=N, frames=frames))


info.sort(key=lambda x: (x['name'], x['celltype'], x['ff'], x['N'], x['frames']))

namelen = max([len(item['name']) for item in info])
cellen = max([len(item['celltype']) for item in info])
Nlen = max([len(item['N']) for item in info])
fflen = max([len(item['ff']) for item in info])

print("Setup frame only:")
print("%*s    %*s    %*s    %*s" %(namelen, "Name", cellen, "celltype", fflen, "ff", Nlen, "N"))
print("-------------------------------------------------")
for f in info:
  if f['frames'] == 0:
    print("%*s    %*s    %*s    %*s" %(namelen, f['name'], cellen, f['celltype'], fflen, f['ff'], Nlen, f['N']))

print("\nAnimations:")
print("%*s    %*s    %*s    %*s    %s" %(namelen, "Name", cellen, "celltype", fflen, "ff", Nlen, "N", "frames"))
print("-----------------------------------------------------------")
for f in info:
  if f['frames'] > 0:
    print("%*s    %*s    %*s    %*s    %6s" %(namelen, f['name'], cellen, f['celltype'], fflen, f['ff'], Nlen, f['N'], f['frames']))
