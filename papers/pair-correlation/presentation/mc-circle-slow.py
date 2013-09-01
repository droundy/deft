from __future__ import division
import time, sys, pylab
from pylab import *
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.animation as animation

if len(sys.argv) > 1 and sys.argv[1] == 'show':
  show = True
else:
  show = False

ff = .4

# setup parameters
x0_setup = 1.2
y0_setup = 1.2
dx_setup = 2.8
dy_setup = 2.5


r = 1
area = pi*r**2

# box
lenx = 20
leny = 10
edge = 1

n = int(round(ff*lenx*leny/area))

# density
dx = 0.1
histogram = zeros(lenx/dx)
xcoords = arange(0, lenx, dx) + dx/2

# Movement
scale = 0.5

def move(pos):
  x = 2*rand()-1
  y = 2*rand()-1
  r2 = x*x + y*y
  while (r2 > 1):
    x = 2*rand()-1
    y = 2*rand()-1
    r2 = x*x + y*y
  fac = sqrt(-2*log(r2)/r2)
  newpos = pos + fac*scale*array([x, y])
  while newpos[1] > leny: newpos[1] -= leny
  while newpos[1] < 0: newpos[1] += leny
  return newpos

def periodic_diff(a, b):
  v = b - a
  while v[1] > leny/2: v[1] -= leny
  while v[1] < -leny/2: v[1] += leny
  return v

def distsq(a, b):
  v = periodic_diff(a, b)
  return v[0]*v[0] + v[1]*v[1]

def touch(a, b):
  if (distsq(a, b) < 4*r*r):
    return True
  return False
##########################


circles = zeros((n+1, 2))
circles[n] = array([-100,-100])
colors = zeros(n+1)

# initial setup
x_setup = x0_setup
y_setup = y0_setup

for i in xrange(n):
  circles[i, 0] = x_setup
  circles[i, 1] = y_setup
  x_setup += dx_setup
  if (x_setup + r > lenx):
    x_setup = x0_setup
    y_setup += dy_setup


# Plotting
fig = figure()
ax = fig.add_subplot(111)

ax.set_xlim(-edge, lenx+edge)
ax.set_ylim(-edge, leny+edge)


#line, = ax2.plot(xcoords, histogram)
#ax2.set_xlim(-edge, lenx+edge)
line2 = ax.arrow(-10, -10, -10, -10)

def init():
  return 0
  #coll = PolyCollection(coords)
  #ax.add_collection(coll)
  #line.set_ydata(histogram)



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

def animate(p):
  global count, circles, skip, colors, histogram, i, line2, success
  if (i >= n):
    i = 0
  keep = True
  temp = move(circles[i])
  if temp[0] - r < 0 or temp[0] + r > lenx:
    keep = False
  else:
    for j in xrange(n):
      if j != i and touch(temp, circles[j]):
        keep = False
        break
  count += 1
  if keep:
    success += 1
    colors[n] = tempc
  else:
    colors[n] = badc
  circles[n] = temp
  colors[i] = highlightc


  ax.cla()
  for j in xrange(n+1):
    for shift in [0, -leny, leny]:
      circ = Circle(circles[j]-(0, shift), r, color=colors[j])
      if shift != 0:
        circ.set_alpha(0.2)
      else:
        circ.set_alpha(0.75)
      ax.add_artist(circ)

  ax.set_title("Attempted moves: %i, Successful moves: %i" %(count, success))
  ax.axhline(y=0, linestyle='--', linewidth=3, color='b')
  ax.axhline(y=leny, linestyle='--', linewidth=3, color='orange')
  ax.axvline(x=0, linestyle='-', color='k', linewidth=3)
  ax.axvline(x=lenx, linestyle='-', color='k', linewidth=3)

  arlen = .5
  arwidth = .5
  for shift in [0, -leny, leny]:
    delta = periodic_diff(circles[i], temp)
    deltax = temp[0]-circles[i,0]
    deltay = temp[1]-circles[i,1]
    if shift == 0:
      alpha = 1
    else:
      alpha = 0.5
    ax.arrow(circles[i,0], circles[i,1]+shift, delta[0],
                     delta[1], head_width=arwidth,
                     head_length=arlen, linewidth=2, facecolor='gray', alpha = alpha)
  fig.tight_layout()
  if keep:
    circles[i] = temp
    colors[i] = goodc
  else:
    colors[i] = rejectc
  if (i == n-1):
    colors = [defaultc]*(n+1)
  i += 1
  if show:
    time.sleep(.3)

fig.tight_layout()
ax.set_aspect('equal')
ax.set_ylim(-2*edge, leny+2*edge)

if show:
  anim = animation.FuncAnimation(fig, animate, init_func=init)
  pylab.show()

else:
  for p in xrange(5*n):
    animate(p)
  count = 0
  success = 0
  print('Saving animation images.')
  for p in xrange(101):
    animate(p)
    savefig("presentation/anim/mc-slow-%03i.pdf" %p)

