#!/usr/bin/python
import glob

names = glob.glob("figs/mc/vertices/*--1.dat")
bases = [name[:-6] for name in names]

info = []
for base in bases:
  frames = len(glob.glob(base+"*")) - 1
  name = base.split('-')
  celltype = name[0].split('/')[-1]
  ff = name[1]
  polyhedron = name[3]
  N = name[4]
  info.append([celltype, polyhedron, ff, N, frames])


info.sort(key=lambda x: (x[1], x[0], x[2], x[3], x[4]))

namelen = max([len(item[1]) for item in info])
cellen = max([len(item[0]) for item in info])
Nlen = max([len(item[3]) for item in info])
fflen = max([len(item[2]) for item in info])

print("Setup frame only:")
print("%*s    %*s    %*s    %*s" %(namelen, "Name", cellen, "celltype", fflen, "ff", Nlen, "N"))
print("------------------------------------------------------------")
for f in info:
  if f[4] == 0:
    print("%*s    %*s    %*s    %*s" %(namelen, f[1], cellen, f[0], fflen, f[2], Nlen, f[3]))

print("\nAnimations:")
print("%*s    %*s    %*s    %*s    %s" %(namelen, "Name", cellen, "celltype", fflen, "ff", Nlen, "N", "frames"))
print("------------------------------------------------------------")
for f in info:
  if f[4] > 0:
    print("%*s    %*s    %*s    %*s    %s" %(namelen, f[1], cellen, f[0], fflen, f[2], Nlen, f[3], f[4]))
