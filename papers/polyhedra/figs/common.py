import os, glob
from numpy import *

# if N=0 is supplied to these functions, then the first appropriate file is used
def get_N(basename):
  # if "vertices" in basename:
  names = glob.glob('%s*' %(basename))
  if len(names) == 0:
    print('No files of the form "%s-*.dat" found.' %(basename))
    return 0
  Ns = []
  for name in names:
    Ns += [int((name.split('-')[4]).split('.')[0])]
  Ns = list(set(Ns))
  print "Using N = %i. Other possible values: " %(Ns[0]),
  for N in Ns[1:]:
    print N,
  print ''
  return Ns[0]

def read_mc_density(ff, poly, N, celltype):
  fname = "figs/mc/%s-%4.2f-density-%s-%i.dat" %(celltype, ff, poly, N)
  print "using", fname
  if (not os.path.isfile(fname)):
    print("\n%s is not a file.\n\nPerhaps you have the wrong number of %ss?" %(fname, poly))
    bad = array([0,10])
    return 0,0
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

def read_mc_pressure(ff, poly, N, celltype):
  fname = "figs/mc/%s-%4.2f-pressure-%s-%i.dat" %(celltype, ff, poly, N)
  print "using", fname
  if (not os.path.isfile(fname)):
    print("\n%s is not a file.\n\nPerhaps you have the wrong number of %ss?" %(fname, poly))
    bad = array([0,10])
    return 0,0
  data = loadtxt(fname)
  dV = data[0,2]
  return dV, data[1:,0], data[1:,1], data[1:,2]

def read_mc_dimensions(ff, poly, N, celltype):
  e, garbage= read_mc_density(ff, poly, N, celltype)
  if e == 0:
    print("getting dimensions from density file didn't work, making some up")
    return ones(3)*5
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
  data = genfromtxt(fname, skip_header=2)
  f = open(fname)
  line = f.readline()
  f.close
  if "iteration" in line:
    iteration =  int(line.split(":")[-1])
  else: iteration = -1
  dim = loadtxt(fname)[0]
  center = data[:, :3]
  verts = data[:, 3:]
  return dim, center, reshape(verts, (N, len(verts[0])/3, 3)), iteration
