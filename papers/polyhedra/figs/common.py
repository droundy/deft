import os
from numpy import *

def read_mc_density(ff, poly, N, celltype='walls'):
  fname = "figs/mc/%s-%4.2f-density-%s-%i.dat" %(celltype, ff, poly, N)
  print "using", fname
  if (not os.path.isfile(fname)):
    print("\n%s is not a file.\n\nPerhaps you have the wrong number of %ss?" %(fname, poly))
    exit(1)
  data = loadtxt(fname)
  xdens = data[:,1]
  ydens = data[:,2]
  zdens = data[:,3]
  x = data[:,0][xdens >= 0]
  y = data[:,0][ydens >= 0]
  z = data[:,0][zdens >= 0]
  xdens = xdens[xdens >= 0]
  ydens = ydens[ydens >= 0]
  zdens = zdens[zdens >= 0]
  return [x, y ,z], [xdens, ydens, zdens]

def read_mc_dimensions(ff, poly, N, celltype='walls'):
  e, garbage= read_mc_density(ff, poly, N, celltype)

  length = empty(3)
  for i in xrange(3):
    de = e[i][2] - e[i][1]
    length[i] = len(e[i])*de
  return length

def check_vertices(ff, shape, N, celltype, f):
  fname = "figs/mc/vertices/%s-%04.2f-vertices-%s-%i-%i.dat" %(celltype, ff, shape, N, f)
  return os.path.isfile(fname)

def read_vertices(ff, shape, N, celltype, f):
  fname = "figs/mc/vertices/%s-%04.2f-vertices-%s-%i-%i.dat" %(celltype, ff, shape, N, f)
  data = loadtxt(fname)
  center = data[:, :3]
  verts = data[:, 3:]
  return center, reshape(verts, (N, len(verts[0])/3, 3))
