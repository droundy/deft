#!/usr/bin/python2
from __future__ import division
import sys
sys.path.insert(0, '../../papers/polyhedra/figs/')

from pylab import *
import mayavi.mlab as mlab
from tvtk.api import tvtk

import animate_polyhedra as anim

# ------------------------------------------------------------------------------
# Generate talk images and then exit
# ------------------------------------------------------------------------------

if __name__ == '__main__':
  print("Generating figures for polyhedra talk")
  mlab.options.offscreen = True
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
    src, shapes = anim.create_plot_data(shape_list, N, True)
    #print shapes

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
  src, shapes = anim.create_plot_data(shape_list, N, True)

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
  src1, shapes1 = anim.create_plot_data(chunk1, 5)
  src2, shapes2 = anim.create_plot_data(chunk2, 5)

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
  # if not args.hide:
  #   surf1, wires1 = plot_tet(src1, .5)
  #   surf2, wires2 = plot_tet(src2, 0)
  #   mlab.view(azimuth=-111, elevation=13, roll=-77, distance=16, focalpoint=(.3, .3, .3))
  #   a = anim()
  #   mlab.show()


