#!/usr/bin/python
from __future__ import division
import scipy, sys, time, os, argparse
from pylab import *

import common


parser = argparse.ArgumentParser(description='Plot density of polyhedra.')
parser.add_argument('ff', metavar='ff', type=float, help='filling fraction')
parser.add_argument('N', metavar='N', type=int, help='number of polyhedra')
parser.add_argument('shape', metavar='shape', help='type of polyhedra',
                    default='cube', choices=['cube', 'tetrahedron', 'truncated_tetrahedron'])
parser.add_argument('-p', '--periodic',help='will use periodic cell - defaults to walls',
                    action='store_true')
parser.add_argument('-s', '--save', help='will save animation frames as pdfs',
                    action='store_true')
parser.add_argument('-f', '--frames', metavar='frames', type=int,
                    help='number of frames to animate', default=10000)
parser.add_argument('--hide', help='will not call show', action='store_true')
parser.add_argument('-o', '--one_frame', metavar='fnumber', type=int, help='Instead of animating, view the specific frame', default=-100)
parser.add_argument('-m', '--mixture', action='store_true',
                    help='signals that the system may have multiple different shapes \
The default assumption is that it does not, which allows faster runtime')
parser.add_argument('--talk', action='store_true', help='instead of performing normal \
function, generate figures for talk')

args = parser.parse_args()

N = args.N
ff = args.ff
polyhedron = args.shape
save = args.save
frames = args.frames

if args.periodic:
  celltype = 'periodic'
else:
  celltype = 'walls'
print celltype

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
  #print("we have %i vertices, %i faces, and the first face has %i vertices." %(len(verts), len(faces), len(faces[0])))
  return array(faces)

import mayavi.mlab as mlab
from tvtk.api import tvtk

figure = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0), size=(1024,768))

if args.talk:
  print("Generating figures for polyhedra talk")
  # input: "../../talks/polyhedra/dat/background.dat" %()
  data = [np.array(map(double, line.split())) for line in
          open("../../talks/polyhedra/dat/background.dat")]
  N = len(data)
  colors = rand(N, 3)
  colors[0] = [1,0,0]
  colors[1] = [0,1,0]
  colors[2] = [0,0,1]
  colors[3] = [1,0,1]
  colors[4] = [.4, 0, .2]
  shape_list = range(N)
  nvertices = sum([(len(data[i])-3)/3 for i in xrange(N)])
  for i in xrange(N):
    data[i] = data[i][3:]
    shape_list[i] = reshape(data[i], (len(data[i])/3, 3))

  for shape, color in zip(shape_list, colors):
    if (max(color) > 1):
      color /= 255
    faces = get_faces(shape)
    mesh = tvtk.PolyData(points=shape, polys=faces)
    src = mlab.pipeline.add_dataset(mesh)
    mlab.pipeline.surface(src, color=tuple(color))
  mlab.view(azimuth=0, elevation=0, distance=20, focalpoint=(0,0,0))
  figure.scene.save("../../talks/polyhedra/figs/background.png")
  if not args.hide:
    mlab.show()
  exit(0)


dim = common.read_mc_dimensions(ff, polyhedron, N, celltype)
print dim

# Initial setup
if args.one_frame > -100:
  f = args.one_frame
else:
  f = 0
centers, shape_array = common.read_vertices(ff, polyhedron, N, celltype, f)

# Get all the data into one Polydata, that is defined by two large arrays
# The first is an N x 3 array of vertex coordines
# The second is an fN x v array where v is the number of vertices per face
# and f is the number of faces per object.
# If the shapes are all the same, then the face array is just repeated copies
# of the face array for one shape, with each copy incremented by the number of
# vertices added.
nvertices = sum([len(shape_array[i]) for i in xrange(N)])
shapes = shape_array.reshape((nvertices, 3))

cvals = linspace(0, 1, N)
# cvals = ones(N)
# cvals[0] = 0
colors = shape_array[0]*0 + cvals[0]
face0 = get_faces(shape_array[0])
faces = face0
vertices_so_far = len(shape_array[0])
for i in xrange(1,N):
  if args.mixture:
    newface = get_faces(shape_array[i]) + vertices_so_far
  else:
    newface = face0 + vertices_so_far
  vertices_so_far += len(shape_array[i])
  faces = vstack((faces, newface))
  colors = vstack((colors, shape_array[i]*0 + cvals[i]))

mesh = tvtk.PolyData(points=shapes, polys=faces)
mesh.point_data.scalars = colors
mesh.point_data.scalars.name = 'colors'
src = mlab.pipeline.add_dataset(mesh)
s = mlab.pipeline.surface(src, colormap='gist_rainbow')

#colormaps: prism, flag, gist_rainbow, gist_ncar



# Draw the cell
mlab.outline(extent=[0, dim[0], 0, dim[1], 0, dim[2]], color=(0,0,0), line_width=3)
if(celltype == 'walls'):
  sheet_points = array([[0, 0, 0], [dim[0], 0, 0], [dim[0], dim[1], 0], [0, dim[1], 0],
                 [0, 0, dim[2]], [dim[0], 0, dim[2]], [dim[0], dim[1], dim[2]],
                 [0, dim[1], dim[2]]])
  sheet_connections = array([[0, 1, 2, 3], [4, 5, 6, 7]])
  sheetmesh = tvtk.PolyData(points=sheet_points, polys=sheet_connections)
  mlab.pipeline.surface(sheetmesh, opacity=.4, color=(.44,.5,.56))

# if(celltype == 'walls'):
#   nbars = 10
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

#   print bar_points
#   print bar_connections
#   bar_src = mlab.pipeline.scalar_scatter(x, y, z, s)
#   bar_src.mlab_source.dataset.lines = bar_connections
#   bars = mlab.pipeline.stripper(bar_src)
#   mlab.pipeline.surface(bars, color=(0,0,0), line_width=2, opacity=.7)

@mlab.animate(delay=30, ui=False)
def anim():
  global f
  while f < frames:
    print f
    if not common.check_vertices(ff, polyhedron, N, celltype, f):
      if save:
        print("All out of dat files.")
        exit(0)
      f = 0
      print("Looping!")
    print f
    shapes[:] = (common.read_vertices(ff, polyhedron, N, celltype, f)[1]).reshape((nvertices, 3))
    src.update()
    yield
    if save:
      figure.scene.save("figs/anim/%s-%4.2f-%s-%i-%i.png" %(celltype, ff, polyhedron, N, f))
    f += 1

if args.one_frame == -100:
  print args.one_frame, f, frames
  a = anim()

if not args.hide:
  mlab.show()
