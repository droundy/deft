#!/usr/bin/python
from __future__ import division
import scipy, sys, time, os, argparse, numpy.random
from pylab import *

import read, arguments
# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------

parser = argparse.ArgumentParser(
  description='Plot density of polyhedra.', parents = [arguments.parser])
parser.add_argument(
  '-d','--delay', metavar='N', type=int, default=10,
  help='time in ms to pause between frames. Minimum 10')

parser.add_argument(
  '-f','--frames', metavar='N', type=int, default=10000,
  help='number of frames to animate')

parser.add_argument(
  '-o','--one_frame', metavar='N', type=int, default=-100,
  help='Instead of animating, view the specific frame')

parser.add_argument(
  '-b','--begin', metavar='N', type=int, default=0,
  help='start at this frame instead of frame 0')

parser.add_argument(
  '-n','--show_only', metavar='F', type=double, default=0,
  help='will only render a fraction of the polyhedra')

parser.add_argument(
  '-k','--skip', metavar='N', type=int, default=1,
  help='only play every N frames')

parser.add_argument(
  '-a','--alpha', metavar='F', type=double,
  help='sets the alpha of the polyhedra', default=1)

parser.add_argument(
  '-v','--save', action='store_true',
  help='will save animation frames as pdfs')

parser.add_argument(
  '-m','--mixture', action='store_true',
  help="""signals that the system may have multiple different shapes.
  The default assumption is that it does not, which allows faster runtime""")

parser.add_argument(
  '-t','--talk', action='store_true',
  help='instead of performing normal function, generate figures for talk')

parser.add_argument('-w','--hide_walls', action='store_true', help='won\'t draw walls')
parser.add_argument('--xmin', type=double, default=0, help='minimum x-value')
parser.add_argument('--xmax', type=double, default=1000, help='maximum x-value')
parser.add_argument('--ymin', type=double, default=0, help='minimum y-value')
parser.add_argument('--ymax', type=double, default=1000, help='maximum y-value')
parser.add_argument('--zmin', type=double, default=0, help='minimum z-value')
parser.add_argument('--zmax', type=double, default=1000, help='maximum z-value')
parser.add_argument(
  '--shiftx', type=double, default=0, help='shift everything over by this many x')
parser.add_argument(
  '--shifty', type=double, default=0, help='shift everything over by this many y')
parser.add_argument(
  '--shiftz', type=double, default=0, help='shift everything over by this many z')
parser.add_argument(
  '--notext', action='store_true', help='won\'t write any text on plot')

args = parser.parse_args()

ff = args.ff
if args.ratio != 0:
  polyhedron = args.shape + "_%05.2f" %args.ratio
else:
  polyhedron = args.shape
save = args.save
frames = args.frames

if args.periodic:
  celltype = 'periodic'
else:
  celltype = 'walls'
print ("Using %s with %s" %(polyhedron, celltype))

if args.N == 0 and not args.talk:
  N = read.get_N("figs/mc/vertices/%s-%4.2f-vertices-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N

if (args.xmin == 0 and args.xmax == 1000 and args.ymin == 0 and args.ymax == 1000 and
   args.zmin == 0 and args.zmax == 1000):
  partial_cell = False
else: partial_cell = True

# ------------------------------------------------------------------------------
# Define functions
# ------------------------------------------------------------------------------
# returns the distance from the point p to the plane defined by 3 points
def dist_to_plane(point, plane):
  p = point-plane[2]
  e1 = plane[0]-plane[2]
  e2 = plane[1]-plane[2]
  n = cross(e1, e2)
  nhat = n/norm(n)
  return dot(nhat, p)

# takes a 2d array and treats it as a list of vertices
# returns a 2d array, that is a list of indices of vertices per face.
# It will pad each list of vertices so that they are a minimum of length pad
# which is useful for mixtures

# new method, look at all possible triplets of vertices, then reject those that have
# other vertices on both sides of the plane that they span
def get_faces(verts, pad = 0, acc = 0.1):
  # find the distance between neighbors. Assumes all neighbors are equidistant
  faces = []
  for i in xrange(len(verts)):
    u = verts[i]
    for j in xrange(i+1, len(verts)):
      v = verts[j]
      for k in xrange(j+1, len(verts)):
        w = verts[k]
        # now make sure we don't have a duplicate
        keep = True
        for face in faces:
          if (i in face) and (j in face) and (k in face):
            keep = False
            break
        if keep:
          plane = vstack((u, v, w))
          has_neg = False
          has_pos = False
          for l in xrange(len(verts)):
            if l != i and l != j and l != k:
              dist = dist_to_plane(verts[l], plane)
              if (dist > acc): has_pos = True
              elif (dist < -acc): has_neg = True
          if (not has_neg) or (not has_pos):
            # this plane is good!
            face = empty(0)
            for l in xrange(len(verts)):
              if abs(dist_to_plane(verts[l], plane)) < acc:
                face = append(face, l)
            faces.append(face)
  # okay we have our faces, but we need to sort them so they'll connect properly
  sfaces = []
  for face in faces:
    sface = array([face[0]])
    for i in xrange(len(face)-1):
      last = sface[-1]
      dmin = 10000
      for j in face:
        if not j in sface:
          dist = norm(verts[last] - verts[j])
          if dist < dmin:
            dmin = dist
            next_neighbor = j
      sface = append(sface, next_neighbor)
    sfaces.append(sface)
  faces = sfaces
  #print("we have %i vertices, %i faces, and the first face has %i vertices." %(len(verts), len(faces), len(faces[0])))
  # enforce that all faces have the same number of points so it can be a
  # 2d array:
  n = max([len(face) for face in faces])
  n = max(n, pad)
  for i in xrange(len(faces)):
    if len(face) < n:
      faces[i] = hstack((faces[i], ones(n-len(faces[i]))*faces[i][-1]))
  return array(faces)


# the neighbor determination is not very robust, and assumes a regular
# polyhedron at the moment
# def get_faces(verts, pad = 0, acc = 0.1):
#   # find the distance between neighbors. Assumes all neighbors are equidistant
#   u = verts[0]
#   dmin = norm(u-verts[1])
#   for j in xrange(2, len(verts)):
#     v = verts[j]
#     d = norm(u - v)
#     if d < dmin: dmin = d
#   faces = []
#   for i in xrange(len(verts)):
#     u = verts[i]
#     for j in xrange(i+1, len(verts)):
#       v = verts[j]
#       if abs(norm(v - u) - dmin) < acc: # they are neighbors
#         for k in xrange(j+1, len(verts)):
#           w = verts[k]
#           if abs(norm(w - v) - dmin) < acc: # we now have two neighbors
#             face_span = vstack((u, v, w))
#             face = array((i, j, k))
#             keep = True
#             for l in xrange(0, len(verts)):
#               if l != i and l != j and l != k:
#                 x = verts[l]
#                 if abs(dist_to_plane(x, face_span)) < acc:
#                   if l < j: keep = False
#                   size = len(face)
#                   face = resize(face, size+1)
#                   face[size] = l
#             if keep:
#               faces.append(face)
#   # okay we have our faces, but we need to sort them so they'll connect properly
#   sfaces = []
#   for face in faces:
#     sface = array([face[0]])
#     for i in xrange(len(face)-1):
#       for j in face:
#         if abs(norm(verts[sface[-1]]-verts[j])-dmin) < acc:
#           if (j in sface) == False:
#             sface = append(sface, j)
#             break
#     sfaces.append(sface)
#   faces = sfaces
#   #print("we have %i vertices, %i faces, and the first face has %i vertices." %(len(verts), len(faces), len(faces[0])))
#   # enforce that all faces have the same number of points so it can be a
#   # 2d array:
#   n = max([len(face) for face in faces])
#   n = max(n, pad)
#   for i in xrange(len(faces)):
#     if len(face) < n:
#       faces[i] = hstack((faces[i], ones(n-len(faces[i]))*faces[i][-1]))
#   return array(faces)

# Get all the data into one Polydata, that is defined by two large arrays
# The first is an N x 3 array of vertex coordines
# The second is an fN x v array where v is the number of vertices per face
# and f is the number of faces per object.
# If the shapes are all the same, then the face array is just repeated copies
# of the face array for one shape, with each copy incremented by the number of
# vertices added.
def create_plot_data(shape_array, ncolors, mixture=False):
  nvertices = sum([len(shape_array[i]) for i in xrange(len(shape_array))])
  shapes = vstack(shape_array)
  N = len(shape_array)
  shapedim = shape_array[0].shape

  cvals = tile(linspace(0, 1, ncolors), ceil(N/ncolors))[:N]
  numpy.random.shuffle(cvals)

  most_vertices = 0
  if mixture:
    for shape in shape_array:
      if len(shape) > most_vertices:
        most_vertices = len(shape)

  face0 = get_faces(shape_array[0], most_vertices)
  faces = empty((0, len(face0[0])))
  colors = empty((0, 3))
  vertices_so_far = 0

  for i in xrange(N):
    if mixture:
      newface = get_faces(shape_array[i], most_vertices) + vertices_so_far
    else:
      newface = face0 + vertices_so_far
    vertices_so_far += len(shape_array[i])

    faces = vstack((faces, newface))
    newcolor = ones(shapedim)*cvals[i] + (2*numpy.random.random_sample(shapedim)-1)*0.1
    colors = vstack((colors, newcolor))

  mesh = tvtk.PolyData(points=shapes, polys=faces)
  mesh.point_data.scalars = colors
  mesh.point_data.scalars.name = 'colors'
  src = mlab.pipeline.add_dataset(mesh)
  return src, shapes

import mayavi.mlab as mlab
from tvtk.api import tvtk

# ------------------------------------------------------------------------------
# Generating talk images and then exit if run with --talk flag
# ------------------------------------------------------------------------------

if args.talk:
  print("Generating figures for polyhedra talk")
  if args.hide:
    # little cubes
    for j in xrange(3):
      f = open("../../talks/polyhedra/dat/tet-%i.dat" %j)
      lines = f.readlines()[2:]
      f.close()
      data = [np.array(map(double, line.split())) for line in lines]
      N = len(data)
      shape_list = range(N)
      nvertices = sum([(len(data[i])-3)/3 for i in xrange(N)])
      for i in xrange(N):
        data[i] = data[i][3:]
        shape_list[i] = reshape(data[i], (len(data[i])/3, 3))
      figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))
      src, shapes = create_plot_data(shape_list, N, True)
      print shapes

      mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1)
      mlab.pipeline.surface(src, representation='wireframe', color = (0, 0, 0))
      mlab.view(azimuth=0, elevation=0, distance=5, focalpoint=(0,0,0))
      figure.scene.save("../../talks/polyhedra/figs/tet-%i.png" %j)
      mlab.clf()
    # Generate background figure
    f = open("../../talks/polyhedra/dat/background.dat", 'r')
    lines = f.readlines()[2:]
    f.close()
    data = [np.array(map(double, line.split())) for line in lines]
    N = len(data)
    shape_list = range(N)
    nvertices = sum([(len(data[i])-3)/3 for i in xrange(N)])
    for i in xrange(N):
      data[i] = data[i][3:]
      shape_list[i] = reshape(data[i], (len(data[i])/3, 3))

    figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))
    src, shapes = create_plot_data(shape_list, N, True)

    mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1)
    mlab.pipeline.surface(src, representation='wireframe', color = (0, 0, 0))

    mlab.view(azimuth=0, elevation=0, distance=20, focalpoint=(0,0,0))
    figure.scene.save("../../talks/polyhedra/figs/background.png")
    mlab.clf()

  # Generate ice structure image

  # parameters to set:
  rO = .6
  dH = .35
  rH = .3
  frames = 4


  figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))
  data = genfromtxt("../../talks/polyhedra/dat/ice-structure.dat", skip_header=2)
  centers = data[:, :3]
  verts = data[:, 3:]
  N = len(centers)
  shape_array = reshape(verts, (N, len(verts[0])/3, 3))


  chunk1 = zeros_like(shape_array[:int(N/2)])
  chunk2 = zeros_like(shape_array[:int(N/2)])

  for i in xrange(int(N/2)):
    chunk1[i] = shape_array[2*i+1]
    chunk2[i] = shape_array[2*i]

  N1 = len(chunk1)
  N2 = len(chunk2)
  src1, shapes1 = create_plot_data(chunk1, 5)
  src2, shapes2 = create_plot_data(chunk2, 5)

  # Now let's try to draw spheres for ices
  oxygens = centers
  Nmolecules = len(oxygens)

  acc = .2
  # we want to map each H20 molecule to its neighbors
  neighbor_distance = 100000
  for center in centers[1:]:
    neighbor_distance = min(neighbor_distance, linalg.norm(centers[0] - center))

  # Plot bonds:
  have_vectors = False
  bond_vec = empty((4,3))
  lines = empty((0,2))
  for i in xrange(N):
    k = 0
    for j in xrange(N):
      #print abs(linalg.norm(centers[j]-centers[i]))
      if abs(linalg.norm(centers[j]-centers[i]) - neighbor_distance) < acc:
        lines = vstack((lines, array([i, j])))
        k += 1
        if (k == 4):
          if not have_vectors:
            # we have a molecule with 4 neighbors, so we can map the vectors
            # among them. Used for placing extra hydrogens that don't have bonds
            # on the edges at the end
            neighbor_bonds = lines[-1:-5:-1]
            for l in xrange(4):
              bond_vec[l] = centers[neighbor_bonds[l,1]] - centers[neighbor_bonds[l,0]]
            have_vectors = True
            break

  # sort and remove duplicates:
  lines = sort(lines, axis=1)
  ncols = lines.shape[1]
  dtype = lines.dtype.descr * ncols
  struct = lines.view(dtype)

  uniq = unique(struct)
  uniq = uniq.view(lines.dtype).reshape(-1, ncols)
  bonds = uniq

  bondcolors = ones_like(centers)
  bondmesh = tvtk.PolyData(points=centers, lines=bonds)
  bondmesh.point_data.scalars = bondcolors
  bondmesh.point_data.scalars.name = 'colors'

  bondsrc = mlab.pipeline.add_dataset(bondmesh)

  # Now to add hydrogen. Each bond gets one hydrogen,
  # and each molecule gets two hydrogens, so let's go from there
  # This is working with with a structured grid of centers,
  # but it could paint itself into a corner with, say, a random arrangement
  bond_h = zeros(len(bonds))
  molecule_h = zeros(len(oxygens))
  hindex = 0
  hydrogens = zeros((len(oxygens)*2, 3))
  while sum(bond_h) < len(bond_h):
    for i in xrange(len(bonds)):
      if bond_h[i] == 0:
        if(molecule_h[bonds[i,0]] < 2):
          i1 = 0
          i2 = 1
        elif(molecule_h[bonds[i,1]] < 2):
          i1 = 1
          i2 = 0
        else:
          print "algorithm broken, finishing anyway"
          i1 = 0
          i2 = 1
        vec = centers[bonds[i,i2]] - centers[bonds[i,i1]]
        uv = vec/linalg.norm(vec)
        loc = centers[bonds[i,i1]] + dH*uv
        hydrogens[hindex] = loc
        hindex += 1
        bond_h[i] += 1
        molecule_h[bonds[i,i1]] += 1
  # now we need to add extra hydrogen to the molecules on the edges
  for i in xrange(len(molecule_h)):
    if molecule_h[i] < 2:
      # first, find a neighbor so we know whether to use bond_vec or -bond_vec:
      for j in xrange(len(bonds)):
        if i in bonds[j]:
          if i == bonds[j,0]:
            neighbor_i = bonds[j,1]
          else:
            neighbor_i = bonds[j,0]
          break
      neighbor_vec = centers[neighbor_i] - centers[i]
      rotate = 1
      for j in xrange(4):
        if linalg.norm(neighbor_vec - bond_vec[j]) < acc:
          break
        elif linalg.norm(neighbor_vec + bond_vec[j]) < acc:
          rotate = -1
          break
      for j in xrange(4):
        # now let's find where our neighbors are, so we can put hydrogens elsewhere
        has_neighbor = False
        for k in xrange(len(centers)):
          if linalg.norm(centers[i] - bond_vec[j] - centers[k]) < acc:
            has_neighbor = True
            break
        if not has_neighbor:
          hydrogens[hindex] = centers[i] + rotate*dH*bond_vec[j]/linalg.norm(bond_vec[j])
          hindex += 1
          molecule_h[i] += 1
          if molecule_h[i] == 2:
            break

  # actual plotting below

  mlab.points3d(oxygens[:,0], oxygens[:,1], oxygens[:,2], color=(.9,1,1), scale_factor=rO, vmin=0, vmax=1, resolution=16)
  mlab.points3d(hydrogens[:,0], hydrogens[:,1], hydrogens[:,2], color=(.3,.3,1), scale_factor=rH, resolution=16)
  mlab.pipeline.surface(bondsrc, line_width=5, color=(1, 0, .2))

  def plot_tet(src, op):
    surf = mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1, opacity=op)
    wireop = op
    if wireop < 1: wireop /= 3
    wires = mlab.pipeline.surface(src, representation='wireframe', color=(0, 0, 0), opacity=wireop)
    return surf, wires

  if args.hide:
    for i in xrange(frames+1):
      op2 = min(i/(frames/2), 1)
      op1 = max(i/(frames/2)-1, 0)
      surf1, wires1 = plot_tet(src1, op1)
      surf2, wires2 = plot_tet(src2, op2)
      mlab.view(azimuth=-149, elevation=12, roll=143, distance=16, focalpoint=(.3, .3, .3))
      #mlab.view(azimuth=-111, elevation=13, roll=-77, distance=16, focalpoint=(.3, .3, .3))
      figure.scene.save("../../talks/polyhedra/figs/ice-structure-%i.png" %i)
      surf1.remove()
      wires1.remove()
      surf2.remove()
      wires2.remove()

  @mlab.animate(delay=1000, ui=False)
  def anim():
    global bondcolors
    while True:
      print mlab.view(), mlab.roll()
      bondcolors[:] = bondcolors*.5
      bondsrc.update()
      yield
  if not args.hide:
    surf1, wires1 = plot_tet(src1, .5)
    surf2, wires2 = plot_tet(src2, 0)
    mlab.view(azimuth=-111, elevation=13, roll=-77, distance=16, focalpoint=(.3, .3, .3))
    a = anim()
    mlab.show()
  exit(0)

# ------------------------------------------------------------------------------
# Inital plot setup
# ------------------------------------------------------------------------------
figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))

if args.one_frame > -100:
  f = args.one_frame
else:
  f = args.begin
dim, centers, shape_array, iteration = read.read_vertices(ff, polyhedron, N, celltype, f)
if partial_cell:
  for i in xrange(N):
    if (centers[i,0] < args.xmin or centers[i,0] > args.xmax or
        centers[i,1] < args.ymin or centers[i,1] > args.ymax or
        centers[i,2] < args.zmin or centers[i,2] > args.zmax):
      shape_array[i] = zeros_like(shape_array[i])
if args.shiftx > 0:
  for i in xrange(len(centers)):
    val = args.shiftx
    centers[i,0] += val
    if centers[i,0] > dim[0]:
      centers[i,0] -= dim[0]
      val -= dim[0]
    shape_array[i,:,0] += val

if args.shifty > 0:
  for i in xrange(len(centers)):
    val = args.shifty
    centers[i,1] += val
    if centers[i,1] > dim[1]:
      centers[i,1] -= dim[1]
      val -= dim[1]
    shape_array[i,:,1] += val

if args.shiftz > 0:
  for i in xrange(len(centers)):
    val = args.shiftz
    centers[i,2] += val
    if centers[i,2] > dim[2]:
      centers[i,2] -= dim[2]
      val -= dim[2]
    shape_array[i,:,2] += val

if partial_cell:
  dim[0] = min(dim[0], args.xmax)
  dim[1] = min(dim[1], args.ymax)
  dim[2] = min(dim[2], args.zmax)
if args.show_only > 0:
  shape_array = shape_array[:N*args.show_only]

src, shapes = create_plot_data(shape_array, 5, args.mixture)
nvertices = len(shapes)
s =  mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1, opacity=args.alpha)

mlab.view(azimuth=30, elevation=65, distance=28, focalpoint=(dim[0]/2,dim[1]/2,dim[2]/2))

words = [polyhedron,
        celltype,
        'ff = %g' %ff,
        'N = %i' %N,
        '                                       ']
text = "\n".join(words)
if not args.notext:
  mlab.text(.8, .9, text, figure=figure, width=.2)
  itertext = mlab.text(.02, .95, "%08i" %iteration, figure=figure, width=.2)

if not args.notext:
  mlab.orientation_axes(figure=figure, name='bob')

# Draw the cell
cell = mlab.outline(extent=[0, dim[0], 0, dim[1], 0, dim[2]], color=(0,0,0), line_width=3)

if(celltype == 'walls' and not args.hide_walls):
  sheet_points = array([[0, 0, 0], [dim[0], 0, 0], [dim[0], dim[1], 0], [0, dim[1], 0],
                 [0, 0, dim[2]], [dim[0], 0, dim[2]], [dim[0], dim[1], dim[2]],
                 [0, dim[1], dim[2]]])
  sheet_connections = array([[0, 1, 2, 3], [4, 5, 6, 7]])
  sheetmesh = tvtk.PolyData(points=sheet_points, polys=sheet_connections)
  mlab.pipeline.surface(sheetmesh, opacity=.6, color=(1,1,1))

# if(celltype == 'walls'and not args.hide_walls):
#   nbars = 11
#   x = tile(repeat(linspace(0, dim[0], nbars), 2), 2)
#   y = tile(array([0, dim[1]]), 2*nbars)
#   z = hstack((zeros(nbars*2), ones(nbars*2)*dim[2]))
#   s = ones(nbars*4)
#   bar_points = zeros((nbars*4, 3))
#   bar_points[:,0] = x
#   bar_points[:,1] = y
#   bar_points[:,2] = z

#   bar_connections = empty((2*nbars, 2))
#   for i in xrange(2*nbars):
#     bar_connections[i,:] = array([2*i, 2*i+1])

#   bar_src = mlab.pipeline.scalar_scatter(x, y, z, s)
#   bar_src.mlab_source.dataset.lines = bar_connections
#   bars = mlab.pipeline.stripper(bar_src)
#   mlab.pipeline.surface(bars, color=(0,0,0), line_width=3, opacity=.7)


# ------------------------------------------------------------------------------
# Animate
# ------------------------------------------------------------------------------

@mlab.animate(delay=args.delay, ui=False)
def anim():
  global f, dim
  while f <= frames:
    if (not read.check_vertices(ff, polyhedron, N, celltype, f)) or f == frames:
      if save:
        print("All done.")
        exit(0)
      f = args.begin
      print("Looping!")
    newdim, newcenters, newshapes, iteration = read.read_vertices(ff, polyhedron, N, celltype, f)
    if not args.notext:
      itertext.set(text="%08i" %iteration)
    if args.show_only>0:
      newshapes = newshapes[:N*args.show_only]
    if partial_cell:
      for i in xrange(N):
        if (newcenters[i,0] < args.xmin or newcenters[i,0] > args.xmax or
            newcenters[i,1] < args.ymin or newcenters[i,1] > args.ymax or
            newcenters[i,2] < args.zmin or newcenters[i,2] > args.zmax):
          newshapes[i] = zeros_like(newshapes[i])

    shapes[:] = newshapes.reshape((nvertices, 3))
    src.update()
    if partial_cell:
      newdim[0] = min(dim[0], args.xmax)
      newdim[1] = min(dim[1], args.ymax)
      newdim[2] = min(dim[2], args.zmax)
    if newdim[0] != dim[0]:
      dim = newdim
      cell = mlab.outline(extent=[0, dim[0], 0, dim[1], 0, dim[2]], color=(0,0,0), line_width=3)
    yield
    if save:
      print("saving figs/anim/%s-%4.2f-%s-%i-%i.png" %(celltype, ff, polyhedron, N, f))
      figure.scene.save("figs/anim/%s-%4.2f-%s-%i-%i.png" %(celltype, ff, polyhedron, N, f))
    f += args.skip

if args.one_frame == -100:
  a = anim()

# if args.notext:
#   mlab.view(azimuth=175, elevation=54, distance=27, focalpoint=(3, 2.6, 11), roll=1.33)
#   figure.scene.save("../../talks/polyhedra/figs/cube-img-%04.2f.png" %ff)

# @mlab.animate(delay=1000, ui=False)
# def anim2():
#   global bondcolors
#   while True:
#     print mlab.view(), mlab.roll()
#     yield
# a = anim2()

if not args.hide:
  mlab.show()
