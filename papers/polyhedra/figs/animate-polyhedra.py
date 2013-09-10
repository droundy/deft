#!/usr/bin/python
from __future__ import division
import scipy, sys, time, os, argparse, numpy.random
from pylab import *

import common


parser = argparse.ArgumentParser(description='Plot density of polyhedra.')
parser.add_argument('ff', metavar='ff', type=float, help='filling fraction')
parser.add_argument('-N', '--N', metavar='N', type=int, default=0,
                    help="""number of polyhedra, if not supplied then the first file
                            with the proper filling fraction will be used""")
parser.add_argument('-s','--shape', metavar='S', default='truncated_tetrahedron',
                    choices=['cube', 'tetrahedron', 'truncated_tetrahedron'],
                    help='type of polyhedra')
parser.add_argument('-d','--delay', metavar='N', type=int, default=10,
                    help='time in ms to pause between frames. Minimum 10')
parser.add_argument('-f','--frames', metavar='N', type=int, default=10000,
                    help='number of frames to animate')
parser.add_argument('-o','--one_frame', metavar='N', type=int, default=-100,
                    help='Instead of animating, view the specific frame')
parser.add_argument('-b','--begin', metavar='N', type=int, default=0,
                    help='start at this frame instead of frame 0')
parser.add_argument('-n','--show_only', metavar='F', type=double, default=0,
                    help='will only render a fraction of the polyhedra')
parser.add_argument('-k','--skip', metavar='N', type=int, default=1,
                    help='only play every N frames')
parser.add_argument('-a','--alpha', metavar='F', type=double,
                    help='sets the alpha of the polyhedra', default=1)

parser.add_argument('-i','--hide', action='store_true', help='will not call show')
parser.add_argument('-p','--periodic', action='store_true',
                    help='will use periodic cell - defaults to walls otherwise')
parser.add_argument('-v','--save', action='store_true',
                    help='will save animation frames as pdfs')
parser.add_argument('-m','--mixture', action='store_true',
                    help="""signals that the system may have multiple different shapes.
                            The default assumption is that it does not, which allows
                            faster runtime""")
parser.add_argument('-t','--talk', action='store_true',
                    help='instead of performing normal function, generate figures for talk')
parser.add_argument('-w','--hide_walls', action='store_true', help='won\'t draw walls')


args = parser.parse_args()

ff = args.ff
polyhedron = args.shape
save = args.save
frames = args.frames

if args.periodic:
  celltype = 'periodic'
else:
  celltype = 'walls'
print ("Using %s with %s" %(polyhedron, celltype))

if args.N == 0 and not args.talk:
  N = common.get_N("figs/mc/vertices/%s-%4.2f-vertices-%s" %(celltype, ff, polyhedron))
  if N == 0:
    exit(1)
else: N = args.N



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

# the neighbor determination is not very robust, and assumes a regular
# polyhedron at the moment
def get_faces(verts, acc = 0.1):
  # find the distance between neighbors. Assumes all neighbors are equidistant
  u = verts[0]
  dmin = norm(u-verts[1])
  for j in xrange(2, len(verts)):
    v = verts[j]
    d = norm(u - v)
    if d < dmin: dmin = d
  faces = []
  for i in xrange(len(verts)):
    u = verts[i]
    for j in xrange(i+1, len(verts)):
      v = verts[j]
      if abs(norm(v - u) - dmin) < acc: # they are neighbors
        for k in xrange(j+1, len(verts)):
          w = verts[k]
          if abs(norm(w - v) - dmin) < acc: # we now have two neighbors
            face_span = vstack((u, v, w))
            face = array((i, j, k))
            keep = True
            for l in xrange(0, len(verts)):
              if l != i and l != j and l != k:
                x = verts[l]
                if abs(dist_to_plane(x, face_span)) < acc:
                  if l < j: keep = False
                  size = len(face)
                  face = resize(face, size+1)
                  face[size] = l
            if keep:
              faces.append(face)
  # okay we have our faces, but we need to sort them so they'll connect properly
  sfaces = []
  for face in faces:
    sface = array([face[0]])
    for i in xrange(len(face)-1):
      for j in face:
        if abs(norm(verts[sface[-1]]-verts[j])-dmin) < acc:
          if (j in sface) == False:
            sface = append(sface, j)
            break
    sfaces.append(sface)
  faces = sfaces
  #print("we have %i vertices, %i faces, and the first face has %i vertices." %(len(verts), len(faces), len(faces[0])))
  # enforce that all faces have the same number of points so it can be a
  # 2d array:
  n = max([len(face) for face in faces])
  for i in xrange(len(faces)):
    if len(face) < n:
      faces[i] = hstack((faces[i], ones(n-len(faces[i]))*faces[i][-1]))
  return array(faces)

import mayavi.mlab as mlab
from tvtk.api import tvtk

figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))
if args.talk:
  print("Generating figures for polyhedra talk")
  # input: "../../talks/polyhedra/dat/background.dat" %()
  f = open("../../talks/polyhedra/dat/background.dat", 'r')
  lines = f.readlines()[2:]
  f.close()
  data = [np.array(map(double, line.split())) for line in lines]
  N = len(data)

  ncols = N
  cvals = tile(linspace(0, 1, ncols), ceil(N/ncols))[:N]

  shape_list = range(N)
  nvertices = sum([(len(data[i])-3)/3 for i in xrange(N)])
  for i in xrange(N):
    data[i] = data[i][3:]
    shape_list[i] = reshape(data[i], (len(data[i])/3, 3))

  for shape, cval in zip(shape_list, cvals):
    shapedim = shape.shape
    color = ones(shapedim)*cval + (2*numpy.random.random_sample(shapedim)-1)*0.15
    faces = get_faces(shape)
    mesh = tvtk.PolyData(points=shape, polys=faces)
    mesh.point_data.scalars = color
    mesh.point_data.scalars.name = 'color'
    src = mlab.pipeline.add_dataset(mesh)
    mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1)
  mlab.view(azimuth=0, elevation=0, distance=20, focalpoint=(0,0,0))
  figure.scene.save("../../talks/polyhedra/figs/background.png")
  if not args.hide:
    mlab.show()
  exit(0)


# Initial setup
if args.one_frame > -100:
  f = args.one_frame
else:
  f = args.begin
dim, centers, shape_array, iteration = common.read_vertices(ff, polyhedron, N, celltype, f)

if args.show_only > 0:
  shape_array = shape_array[:N*args.show_only]


mlab.view(azimuth=30, elevation=65, distance=20, focalpoint=(0,0,0))


# Get all the data into one Polydata, that is defined by two large arrays
# The first is an N x 3 array of vertex coordines
# The second is an fN x v array where v is the number of vertices per face
# and f is the number of faces per object.
# If the shapes are all the same, then the face array is just repeated copies
# of the face array for one shape, with each copy incremented by the number of
# vertices added.
nvertices = sum([len(shape_array[i]) for i in xrange(len(shape_array))])
shapes = shape_array.reshape((nvertices, 3))


shapedim = shape_array[0].shape

ncols = 5
cvals = tile(linspace(0, 1, ncols), ceil(N/ncols))[:N]

face0 = get_faces(shape_array[0])
faces = empty((0, len(face0[0])))
colors = empty((0, 3))
vertices_so_far = 0
for i in xrange(len(shape_array)):
  if args.mixture:
    newface = get_faces(shape_array[i]) + vertices_so_far
  else:
    newface = face0 + vertices_so_far
  vertices_so_far += len(shape_array[i])
  faces = vstack((faces, newface))
  newcolor = ones(shapedim)*cvals[i] + (2*numpy.random.random_sample(shapedim)-1)*0.15
  colors = vstack((colors, newcolor))

mesh = tvtk.PolyData(points=shapes, polys=faces)
mesh.point_data.scalars = colors
mesh.point_data.scalars.name = 'colors'
src = mlab.pipeline.add_dataset(mesh)
s = mlab.pipeline.surface(src, colormap='jet', vmin=0, vmax=1, opacity=args.alpha)

#colormaps: prism, flag, gist_rainbow, gist_ncar
mlab.orientation_axes(figure=figure, name='bob')
words = [polyhedron,
        celltype,
        'ff = %g' %ff,
        'N = %i' %N,
        '                                       ']
text = "\n".join(words)
mlab.text(.8, .9, text, figure=figure, width=.2)
itertext = mlab.text(.02, .95, "%08i" %iteration, figure=figure, width=.2)



# Draw the cell
cell = mlab.outline(extent=[0, dim[0], 0, dim[1], 0, dim[2]], color=(0,0,0), line_width=3)

# if(celltype == 'walls' and not args.hide_walls):
#   sheet_points = array([[0, 0, 0], [dim[0], 0, 0], [dim[0], dim[1], 0], [0, dim[1], 0],
#                  [0, 0, dim[2]], [dim[0], 0, dim[2]], [dim[0], dim[1], dim[2]],
#                  [0, dim[1], dim[2]]])
#   sheet_connections = array([[0, 1, 2, 3], [4, 5, 6, 7]])
#   sheetmesh = tvtk.PolyData(points=sheet_points, polys=sheet_connections)
#   mlab.pipeline.surface(sheetmesh, opacity=.3, color=(.44,.5,.56))

if(celltype == 'walls'and not args.hide_walls):
  nbars = 11
  x = tile(repeat(linspace(0, dim[0], nbars), 2), 2)
  y = tile(array([0, dim[1]]), 2*nbars)
  z = hstack((zeros(nbars*2), ones(nbars*2)*dim[2]))
  s = ones(nbars*4)
  bar_points = zeros((nbars*4, 3))
  bar_points[:,0] = x
  bar_points[:,1] = y
  bar_points[:,2] = z

  bar_connections = empty((2*nbars, 2))
  for i in xrange(2*nbars):
    bar_connections[i,:] = array([2*i, 2*i+1])

  bar_src = mlab.pipeline.scalar_scatter(x, y, z, s)
  bar_src.mlab_source.dataset.lines = bar_connections
  bars = mlab.pipeline.stripper(bar_src)
  mlab.pipeline.surface(bars, color=(0,0,0), line_width=3, opacity=.7)

if args.delay >= 10:
  delay = args.delay
else:
  delay = 10
@mlab.animate(delay=delay, ui=False)
def anim():
  global f, dim
  while f < frames:
    if not common.check_vertices(ff, polyhedron, N, celltype, f):
      if save:
        print("All out of dat files.")
        exit(0)
      f = args.begin
      print("Looping!")
    newdim, newcenters, newshapes, iteration = common.read_vertices(ff, polyhedron, N, celltype, f)
    itertext.set(text="%08i" %iteration)
    if args.show_only>0:
      newshapes = newshapes[:N*args.show_only]
    shapes[:] = newshapes.reshape((nvertices, 3))
    src.update()
    if newdim[0] != dim[0]:
      dim = newdim
      cell = mlab.outline(extent=[0, dim[0], 0, dim[1], 0, dim[2]], color=(0,0,0), line_width=3)
    #cell.mlab_source.extent=[0, dim[0], 0, dim[1], 0, dim[2]]
    #cell.set(xmax=dim[0])
    yield
    if save:
      print("saving figs/anim/%s-%4.2f-%s-%i-%i.png" %(celltype, ff, polyhedron, N, f))
      figure.scene.save("figs/anim/%s-%4.2f-%s-%i-%i.png" %(celltype, ff, polyhedron, N, f))
    f += args.skip

if args.one_frame == -100:
  a = anim()

if not args.hide:
  mlab.show()
