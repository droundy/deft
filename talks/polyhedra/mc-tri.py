#!/usr/bin/python
from __future__ import division
import time, sys, matplotlib, numpy
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.animation as animation

ff = .5

cdefault = .3
cspecial = .9

# setup parameters
x0_setup = 1
y0_setup = 1
dx_setup = 2
dy_setup = 1.1
dy2_setup = 1.2

# for triangles
l = 1 # length from center to vertex
xdiff = 0.5
ydiff = 1/(2*sqrt(3))
ydiff2 = 1/sqrt(3)
area = (sqrt(3)/2*l)*(l/2)/2*6

# box
lenx = 20
leny = 10
edge = 1

n = int(ff*lenx*leny/area)


# density
dx = 0.25
dtheta = 2*pi/3/lenx*dx
histogram = zeros(round(lenx/dx))
angle_histogram = zeros((round(2*pi/3/dtheta), round(lenx/dx)))

xcoords = arange(0, lenx, dx) + dx/2

# Movement
phi_m = pi/32
r_m = .1

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
  ang = angle + fac*phi_m*x
  while ang > 2*pi/3:
    ang -= 2*pi/3
  while ang < 0:
    ang += 2*pi/3
  return ang

# Functions
def make_square(xs, ys):
  dx = xs[2] - xs[1]
  new_ys = repeat(ys, 2)
  new_xs = repeat(xs, 2)
  new_xs += (dx/2.0)
  new_xs[::2] -= dx
  return new_xs, new_ys

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
  for i in xrange(3):
    p1 = v1[i-1]
    p2 = v1[i]
    m1 = (p2[1]-p1[1])/(p2[0]-p1[0])
    b1 = p1[1]-m1*p1[0]
    for j in xrange(3):
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
##########################


centers = zeros((n+1, 2))
centers[n] = array([-100,-100])
angles = zeros(n+1)
colors = zeros(n+1)

# initial setup
x_setup = x0_setup
y_setup = y0_setup
phi_setup = 0
for i in xrange(n):
  centers[i, 0] = x_setup
  centers[i, 1] = y_setup
  angles[i] = phi_setup

  if max(verts(centers[i], angles[i])[:,0]) > lenx:
    ang = pi/3
    x_setup = x0_setup + (ang - phi_setup)/ang*x0_setup
    y_setup += (ang - phi_setup)/ang*dy_setup + phi_setup/ang*dy2_setup
    phi_setup = ang - phi_setup
    centers[i, 0] = x_setup
    centers[i, 1] = y_setup
    angles[i] = phi_setup
  x_setup += dx_setup

coords = zeros((n, 3, 2))
for i in xrange(n):
  coords[i] = verts(centers[i], angles[i])
# Plotting
fig = figure()
ax = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharex=ax)
ax3 = fig.add_subplot(223, sharex=ax)
ax4 = fig.add_subplot(224, sharex=ax, sharey=ax)


ax3.set_ylim(-1/3, 1/3)
ax3.set_yticks((-1/3, -1/6, 0, 1/6, 1/3))
ax3.set_yticklabels((r'$-\pi/3$', r'$-\pi/6$', '0', r'$\pi/6$', r'$\pi/3$'))
ax3.set_yticklabels(['   .   ']*5)

ax5 = ax3.twiny()
ax5.set_xticklabels([])
#ax5.set_yticklabels([])
ax5.set_aspect('equal')
#print ax5.set_ylim(ax3.get_ylim())

# rotated triangle pictures:
pic_x = -.03
pic_l = .04
pic_angles = array([-pi/3, -pi/6, 0, pi/6, pi/3])
pic_centers = array([[pic_x, i] for i in pic_angles/pi])

temp = l
l = pic_l
pic_coords = [verts(pic_centers[i], -pic_angles[i]) for i in xrange(len(pic_centers))]
l = temp

coll = PolyCollection(pic_coords, clip_on=False, color=cm.jet(cdefault))
ax5.add_collection(coll)

# diagrams on fourth plot:
# temp = l
# l = 3
# tri = array((3,0))
# tri_coords = verts(tri, 0)
# ang = -pi/6
# new_tri_coords = verts(tri, ang)
# l = temp

ax4.axis('off')
# ax4.add_patch(Polygon(tri_coords, fill=False))
# ax4.add_patch(Polygon(new_tri_coords, alpha = .5, color=cm.jet(cdefault)))

# line_len = 7
# ax4.plot([tri[0], tri[0]], [tri[1], tri[1]+line_len], '-k')
# ax4.plot([tri[0], tri[0]+line_len*sin(ang)], [tri[1], tri[1]+line_len*cos(ang)], '--k')

# phis = arange(0, ang, -.01)
# circ_r = 4
# ax4.plot(tri[0] + circ_r*sin(phis), tri[1] + circ_r*cos(phis), '-k')

# ax4.annotate('$\\varphi$', xy=(tri[0] + 6*sin(ang/2)-.5, tri[1] + 6*cos(ang/2)), fontsize=20)


ax4.annotate(r'$n(x) = \frac{\mathrm{counts\ in\ slice}}{\mathrm{total\ counts}}\cdot\frac{N}{A}$', xy=(2, 5), fontsize=20)
#fig, (ax, ax2, ax3) = subplots(3,2, sharex=True)
#fig2 = figure(2)
#ax3 = fig2.add_subplot(111)

ax.axvline(x=0, linestyle='-', color='k', linewidth=3)
ax.axvline(x=lenx, linestyle='-', color='k', linewidth=3)
ax.axhline(y=0, linestyle='--', linewidth=3, color='b')
ax.axhline(y=leny, linestyle='--', linewidth=3, color='orange')



ax.set_ylim(-edge, leny+edge)
ax.set_xlim(-edge, lenx+edge)
ax2.set_title("number density, $n$")
ax3.set_title("angle")

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax3.get_xaxis().set_visible(False)
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)


xnew, density = make_square(xcoords, histogram)
x = arange(0, lenx, dx)
theta = arange(-1/3, 1/3, dtheta/pi)
x,theta = meshgrid(x, theta)

line, = ax2.plot(xnew, density)
#angplot = ax3.pcolormesh(x, theta, angle_histogram, vmin=0, vmax=.01)
ax2.set_xlim(-edge, lenx+edge)


def init():
  coll = PolyCollection(coords)
  ax.add_collection(coll)
  line.set_ydata(density)
  return line, ax2

count = 0
success = 0
skip = 150

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

def mc():
  global count, centers, angles, coords, skip, colors, histogram, success, angle_histogram
  for j in xrange(skip):
    count += 1
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
      if keep:
        success += 1
    # add histogram counts:
    bins = linspace(0, lenx, len(histogram)+1)
    histogram += numpy.histogram(centers[:,0], bins=bins)[0]
    xmax = 0
    thetamax = 0
    for i in xrange(n):
      x = int(centers[i,0]/dx)
      theta = int(angles[i]/dtheta)
      xmax = max(x, xmax)
      thetamax = max(theta, thetamax)
      angle_histogram[theta, x] += 1

#fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.517, 0.05, 0.02, 0.4])

def animate(p):
  global count, centers, angles, coords, skip, colors, histogram, success, angle_histogram
  mc()
  density = histogram/leny/dx/count
  xnew, density = make_square(xcoords, density)
  line.set_ydata(density)
  ang_vals = angle_histogram/count/pi
  ax3.collections = []
  angplot = ax3.pcolormesh(x, theta, ang_vals, vmin=0, vmax=.01, edgecolors='face', cmap=cm.hot_r)
  #angplot = ax3.contourf(x, theta, ang_vals, levels=arange(0, .0105, .001), extend="max", rasterized=True)
  cbar_ax.collections = []
  cs = fig.colorbar(angplot, cax=cbar_ax, ticks=[])
  cs.cmap.set_over('k')
  cs.set_clim([0, .01])
  ax2.set_ylim(0, 0.8)
  for i in xrange(n):
    coords[i] = verts(centers[i], angles[i])
  coll = PolyCollection(coords)

  colors = zeros(n) + cdefault
  colors[0] = cspecial
  coll.set_color([cm.jet(val) for val in colors])
  ax.collections=[]
  ax.add_collection(coll)
  ax.set_title("Attempted: %6i, Successful: %6i" %(count*n, success))
  #fig.tight_layout()
  print p
  return line, ax2, ax3

fig.tight_layout()
ax.set_aspect('equal')
ax.set_ylim(-edge, leny+edge)

#fig.set_size_inches(12,7)

if 'show' in sys.argv:
  anim = animation.FuncAnimation(fig, animate, init_func=init)
  show()
  exit(0)

print("Generating mc-tri animation images for polyhedra paper.")
for i in xrange(30):
  initialize()
for p in xrange(30):
  animate(p)
  savefig("anim/mc100-%4.2f-%03i.png" %(ff, p))
