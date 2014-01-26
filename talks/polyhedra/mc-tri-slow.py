#!/usr/bin/python
from __future__ import division
import time, sys, matplotlib
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.animation as animation

ff = .2

# setup parameters
x0_setup = 1
y0_setup = 1
dx_setup = 2.4
dy_setup = 2
dy2_setup = 2.7

# for triangles
l = 1 # length from center to vertex
xdiff = 0.5
ydiff = 1/(2*sqrt(3))
ydiff2 = 1/sqrt(3)
area = (sqrt(3)/2*l)*(l/2)/2*6

# box
lenx = 20
leny = 10
xedge = 1.5
yedge = 6

n = int(ff*lenx*leny/area)


# density
dx = 0.1
histogram = zeros(lenx/dx)
xcoords = arange(0, lenx, dx) + dx/2

# Movement
phi_m = pi/8
r_m = .5

def move(pos):
  x = 2*rand()-1
  y = 2*rand()-1
  r2 = x*x + y*y
  while (r2 > 1):
      x = 2*rand()-1
      y = 2*rand()-1
      r2 = x*x + y*y
  fac = sqrt(-2*log(r2)/r2)
  newpos = pos + fac*r_m*array([x, y])
  while newpos[1] > leny: newpos[1] -= leny
  while newpos[1] < 0: newpos[1] += leny
  return newpos

def turn(angle):
  x = 2*rand()-1
  y = 2*rand()-1
  r2 = x*x + y*y
  while (r2 > 1):
      x = 2*rand()-1
      y = 2*rand()-1
      r2 = x*x + y*y
  fac = sqrt(-2*log(r2)/r2)
  return angle + fac*phi_m*x

# Functions
def distsq(a, b):
  v = periodic_diff(a, b)
  return v[0]*v[0] + v[1]*v[1]

def periodic_diff(a, b):
  v = b - a
  while v[1] > leny/2: v[1] -= leny
  while v[1] < -leny/2: v[1] += leny
  return v


def touch(pos1, pos2, phi1, phi2):
  if (distsq(pos2, pos1) > 4*l*l):
    return False
  if (distsq(pos2, pos1) < l*l):
    return True
  v1 = verts([0,0], phi1)
  v2 = verts(periodic_diff(pos1, pos2), phi2)
  for i in range(3):
    p1 = v1[i-1]
    p2 = v1[i]
    m1 = (p2[1]-p1[1])/(p2[0]-p1[0])
    b1 = p1[1]-m1*p1[0]
    for j in range(3):
      p3 = v2[j-1]
      p4 = v2[j]
      m2 = (p4[1]-p3[1])/(p4[0]-p3[0])
      b2 = p3[1]-m2*p3[0]
      xint = (b1-b2)/(m2-m1)
      if xint < p1[0] and xint > p2[0]:
        if xint < p3[0] and xint > p4[0]:
          return True
        if xint < p4[0] and xint > p3[0]:
          return True
      if xint < p2[0] and xint > p1[0]:
        if xint < p3[0] and xint > p4[0]:
          return True
        if xint < p4[0] and xint > p3[0]:
          return True
  return False

def verts(pos, phi):
  x = pos[0]
  y = pos[1]
  return array([[x+l*sin(phi), y+l*cos(phi)], [x+l*sin(phi+2*pi/3), y+l*cos(phi+2*pi/3)], [x+l*sin(phi + 4*pi/3), y+l*cos(phi + 4*pi/3)]])

def chop(triangle):
  if min(triangle[:,1]) > 0 and max(triangle[:,1]) < leny:
    shape = vstack((triangle, triangle[0]))
    offset = zeros_like(shape)
    offset[:,1] = leny
    return shape, shape, shape+offset, shape-offset

  sub = zeros_like(triangle)
  top = False
  if min(triangle[:,1]) > leny/2:
    sub[:,1] = leny
    top = True

  tri = triangle - sub
  if median(tri[:,1]) < 0:
    loner_index = argmax(tri[:,1])
  else:
    loner_index = argmin(tri[:,1])
  loner = zeros(2)
  duo1 = zeros(2)
  duo2 = zeros(2)
  for i in range(3):
    if i == loner_index:
      loner = tri[i]
    elif duo1[1] == 0:
      duo1 = tri[i]
    else:
      duo2 = tri[i]
  xint1 = loner[0] - (duo1[0]-loner[0])/(duo1[1]-loner[1])*loner[1]
  xint2 = loner[0] - (duo2[0]-loner[0])/(duo2[1]-loner[1])*loner[1]

  p1 = array([xint1, sub[0,1]])
  p2 = array([xint2, sub[0,1]])
  shape1 = vstack((p1, loner+sub[0], p2, p1))
  shape2 = vstack((p1, duo1+sub[0], duo2+sub[0], p2))
  offset = zeros_like(shape1)
  offset[:,1] = leny
  if loner[1] < 0:
    if top:
      return shape1, shape2-offset, shape1-offset, shape2
    else:
      return shape2, shape1+offset, shape1, shape2+offset
  else:
    if top:
      return shape2, shape1-offset, shape2-offset, shape1
    else:
      return shape1, shape2+offset, shape2, shape1+offset
##########################


centers = zeros((n+1, 2))
centers[n] = array([-100,-100])
angles = zeros(n+1)
colors = zeros(n+1)

# initial setup
x_setup = x0_setup
y_setup = y0_setup
phi_setup = 0
for i in range(n):
  centers[i, 0] = x_setup
  centers[i, 1] = y_setup
  angles[i] = phi_setup

  if max(verts(centers[i], angles[i])[:,0]) > lenx:
    x_setup = x0_setup + (pi - phi_setup)/pi*x0_setup
    y_setup += (pi - phi_setup)/pi*dy_setup + phi_setup/pi*dy2_setup
    phi_setup = pi - phi_setup
    centers[i, 0] = x_setup
    centers[i, 1] = y_setup
    angles[i] = phi_setup
  x_setup += dx_setup

coords = zeros((n+1, 4, 2))
coords2 = zeros((n+1, 4, 2))
extras = zeros((n+1, 4, 2))
extras2 = zeros((n+1, 4, 2))

# Plotting
fig = figure()
ax = fig.add_subplot(111)


line2 = ax.arrow(-10, -10, -10, -10)

def init():
  coll = PolyCollection(coords)
  ax.add_collection(coll)
  #line.set_ydata(histogram)
  return ax, line2


alpha = .75
defaultc = (0,.7,1,1)
highlightc = (0,0,1,1)
goodc = (.3,1,.5,1)
badc = (1,0,0,1)
rejectc = (.6,.1,.1,1)
tempc = (.1,1,.1,1)
count = 0
success = 0
skip = 1
i = 0

colors = [defaultc]*(n+1)
def initialize():
  global centers, angles
  for i in xrange(n):
    keep = True
    temp = move(centers[i])
    tempa = turn(angles[i])
    xs = verts(temp, tempa)[:,0]
    if max(xs) > lenx or min(xs) < 0:
      keep = False
    else:
      for j in xrange(n):
        if j != i and touch(temp, centers[j], tempa, angles[j]):
          keep = False
          break
    if keep:
      centers[i] = temp
      angles[i] = tempa

def animate(p):
  global count, centers, angles, coords, extras, skip, colors, histogram, i, line2, success
  if (i >= n):
    i = 0
  keep = True
  ax.cla()
  temp = move(centers[i])
  tempa = turn(angles[i])
  coords[n], coords2[n], extras[n], extras2[n] = chop(verts(temp, tempa))
  xs = verts(temp, tempa)[:,0]
  if max(xs) > lenx or min(xs) < 0:
    keep = False
  else:
    for j in xrange(n):
      if j != i and touch(temp, centers[j], tempa, angles[j]):
        keep = False
        break

  count += 1
  if keep:
    success += 1
    colors[n] = tempc
  else:
    colors[n] = badc


  # for k in range(n+1):
  #   for j in range(len(extras[k,:])):
  #     lines = plot([extras[k,j-1,0], extras[k,j,0]], [extras[k,j-1,1], extras[k,j,1]], '--', color=colors[k], linewidth=2, alpha=alpha)
  #     lines = plot([extras2[k,j-1,0], extras2[k,j,0]], [extras2[k,j-1,1], extras2[k,j,1]], '--', color=colors[k], linewidth=2, alpha=alpha)

  colors[i] = highlightc
  # swap things around so the moving triangles are always on top
  c = concatenate((coords, coords2), axis=0)
  c_extra = concatenate((extras, extras2), axis=0)

  cols = colors + colors

  swap = c[-2]
  swap_extra = c_extra[-2]
  c[-2] = c[n]
  c_extra[-2] = c_extra[n]
  c[n] = swap
  c_extra[n] = swap_extra

  swap = cols[-2]
  cols[-2] = cols[n]
  cols[n] = swap

  coll = PolyCollection(c)

  coll.set_color(cols)
  coll.set_alpha(alpha)

  coll_extra = PolyCollection(c_extra)
  coll_extra.set_color(cols)
  coll_extra.set_alpha(.2)
  #coll.set_color([cm.jet(x) for x in cols])

  ax.collections=[]
  ax.add_collection(coll_extra)
  ax.add_collection(coll)
  ax.set_title("Attempted moves: %i, Successful moves: %i" %(count, success))
  ax.axhline(y=0, linestyle='--', linewidth=3, color='b', zorder=1)
  ax.axhline(y=leny, linestyle='--', linewidth=3, color='orange', zorder=1)
  ax.axvline(x=0, linestyle='-', color='k', linewidth=3, zorder=1)
  ax.axvline(x=lenx, linestyle='-', color='k', linewidth=3, zorder=1)

  arlen = .5
  arwidth = .5
  delta = periodic_diff(centers[i], temp)
  for shift in [0, -leny, leny]:
    if shift == 0:
      aralpha = 1
    else:
      aralpha = .5
    line2 = ax.arrow(centers[i,0], centers[i,1]+shift, delta[0],
                     delta[1], head_width=arwidth, head_length=arlen,
                     linewidth=2, facecolor='slategray', zorder=3, alpha=aralpha)
  fig.tight_layout()
  if keep:
    centers[i] = temp
    angles[i] = tempa
    coords[i], coords2[i], extras[i], extras2[i] = chop(verts(centers[i], angles[i]))
    colors[i] = goodc
  else:
    colors[i] = rejectc
  if (i == n-1):
    colors = [defaultc]*(n+1)
  i += 1
  return ax, line2

fig.tight_layout()
ax.set_aspect('equal')

ax.set_xlim(-xedge, lenx+xedge)
ax.set_ylim(-yedge, leny+yedge)
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)

print("Generating mc-tri-slow animation images for polyhedra paper.")
for p in xrange(50):
  initialize()

count = 0
success = 0
for i in xrange(n+1):
  coords[i], coords2[i], extras[i], extras2[i] = chop(verts(centers[i], angles[i]))

# anim = animation.FuncAnimation(fig, animate, init_func=init)
# fig.set_size_inches(8,4)
# show()

for p in xrange(30):
  animate(p)
  savefig("anim/mc-slow-%03i.pdf" %p)

