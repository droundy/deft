#!/usr/bin/python2
from __future__ import division
import scipy, sys, time, os, argparse, numpy.random
from pylab import *

import read, arguments
import mayavi.mlab as mlab
from tvtk.api import tvtk

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

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
if __name__ == '__main__':
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

  if args.hide:
    mlab.options.offscreen = True

  ff = args.ff
  if args.ratio != 0:
    polyhedron = args.shape + "_%05.2f" %args.ratio
  else:
    polyhedron = args.shape
  frames = args.frames

  if args.periodic:
    celltype = 'periodic'
  else:
    celltype = 'walls'
  print ("Using %s with %s" %(polyhedron, celltype))

  if args.N == 0:
    N = read.get_N("figs/mc/vertices/%s-%4.2f-vertices-%s" %(celltype, ff, polyhedron))
    if N == 0:
      exit(1)
  else: N = args.N

  if (args.xmin == 0 and args.xmax == 1000 and args.ymin == 0 and args.ymax == 1000 and
      args.zmin == 0 and args.zmax == 1000):
    partial_cell = False
  else: partial_cell = True


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
    if iteration == -1:
      iterword = 'initial grid'
    elif iteration == -2:
      iterword = 'post initialization'
    elif iteration == -3:
      iterword = 'latest save'
    else:
      iterword = '%08i' %iteration
    itertext = mlab.text(.02, .95, iterword, figure=figure, width=.2)

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
        if args.save:
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
      # fixme : yield screws things up even when args.hide == True
      # if not args.hide:
        # yield
      if args.save:
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

  mlab.show()
