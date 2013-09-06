#!/usr/bin/python

from __future__ import division
import time, sys, matplotlib

if not ('show' in sys.argv):
  matplotlib.use('Agg')

from pylab import *
seed(3) # seed random number generator
import pylab
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.animation as animation

ff = .4


r = 1
delta_contact = 0.03
area = pi*r**2

# box
lenx = 20
leny = 10
edge = 1

N = int(round(ff*lenx*leny/area))
plotmoves = 1
if 'pair' in sys.argv or 'density' in sys.argv:
  plotmoves = 30

# density
dx = 0.2
histogram = zeros(lenx/dx)

xcoords_doubled = zeros(2*len(histogram))
for i in xrange(len(histogram)):
  xcoords_doubled[2*i] = i*dx
  xcoords_doubled[2*i+1] = (i+1)*dx
histogram_doubled = zeros_like(xcoords_doubled)
contact_histogram = zeros_like(xcoords_doubled)
nA_histogram = zeros_like(xcoords_doubled)
contact_counts = 0

pairdx = 0.5
xxx = arange(0, lenx+pairdx/2, pairdx)
yyy = arange(0, leny+pairdx/2, pairdx)
XXX, YYY = meshgrid(xxx, yyy)
pairdensity = zeros_like(XXX)

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

def almost_touch(a, b):
  if (distsq(a, b) < 4*(r+delta_contact)*(r+delta_contact)):
    return True
  return False
##########################


circles = zeros((N+1, 2))
circles[N] = array([-100,-100])
colors = zeros(N+1)

# setup parameters
dx_setup = sqrt((lenx-2)*leny/N)
nx_setup = ceil((lenx-2)/dx_setup)
ny_setup = ceil(leny/dx_setup)
n_setup = nx_setup*ny_setup
x_setup = linspace(dx_setup/2, lenx-dx_setup/2, nx_setup)
y_setup = linspace(dx_setup/2, leny-dx_setup/2, ny_setup)
X_setup, Y_setup = meshgrid(x_setup, y_setup)

occupied_setup = ones((nx_setup,ny_setup))
while sum(occupied_setup) > N:
  occupied_setup[floor(nx_setup*rand()), floor(ny_setup*rand())] = 0

i_setup = 0
for i in xrange(int(nx_setup)):
  for j in xrange(int(ny_setup)):
    if occupied_setup[i,j]:
      circles[i_setup, 0] = x_setup[i]
      circles[i_setup, 1] = y_setup[j]
      i_setup += 1
if 'pair' in sys.argv:
  oldzerox = circles[0,0]
  oldzeroy = circles[0,1]
  circles[0,0] = 1
  circles[0,1] = leny/2
  for i in xrange(1, N):
    if touch(circles[i], circles[0]):
      circles[i,0] = oldzerox
      circles[i,1] = oldzeroy
      break

# Plotting
fig = figure()
if 'pair' in sys.argv:
  ax = fig.add_subplot(111)
  ax.xaxis.set_visible(False)
else:
  ax2 = fig.add_subplot(111)
  ax2.xaxis.set_visible(False)
  ax = ax2.twinx()
if 'gsigma' in sys.argv:
  ax2.set_yticks([])
ax.yaxis.set_visible(False)

ax.set_xlim(-edge, lenx+edge)
ax.set_ylim(-edge, leny+edge)

if 'density' in sys.argv:
  line, = ax2.plot(xcoords_doubled, histogram_doubled)
  ax2.set_ylabel('density $n(x)$')
  ax2.axvline(x=0, linestyle='-', color='k', linewidth=3)
  ax2.axvline(x=lenx, linestyle='-', color='k', linewidth=3)
  ax2.set_ylim(0, 2)
elif 'gsigma' in sys.argv:
  line, = ax2.plot(xcoords_doubled, contact_histogram)
  ax2.set_ylabel('$g_\sigma$')
  ax2.axvline(x=0, linestyle='-', color='k', linewidth=3)
  ax2.axvline(x=lenx, linestyle='-', color='k', linewidth=3)
  ax2.set_ylim(0, 2)
elif 'pair' in sys.argv:
  pairdensity[0,0] = 0.01
  pcolormesh(XXX, YYY, (pairdensity.max()-pairdensity)/sum(pairdensity), cmap='hot')

def init():
  return 0

def setup():
  global circles, histogram_doubled
  for i in xrange(1,N):
    keep = True
    temp = move(circles[i])
    if temp[0] - r < 0 or temp[0] + r > lenx:
      keep = False
    else:
      for j in xrange(N):
        if j != i and touch(temp, circles[j]):
          keep = False
          break
    if keep: circles[i] = temp
    if 'gsigma' in sys.argv:
      binnum = round(circles[i][0]/dx) % len(histogram)
      histogram_doubled[2*binnum] += 1 # add another count to histogram
      histogram_doubled[2*binnum+1] += 1 # add another count to histogram
      for j in xrange(len(histogram)):
        xjp = j*dx + dx
        xjm = j*dx
        if abs(xjp - circles[i][0]) <= 2 and abs(xjm - circles[i][0]) <= 2:
          lineseg = abs(sqrt(4 - (xjp-circles[i][0])**2) -  sqrt(4 - (xjm-circles[i][0])**2))
          nA_histogram[2*j] += lineseg
          nA_histogram[2*j+1] += lineseg

defaultc = (0,.7,1,1)
highlightc = (0,0,1,1)
contactc = (0,0,0.25)
if 'pair' in sys.argv:
  goodc = defaultc
  rejectc = defaultc
  highlightc = defaultc
else:
  goodc = (.3,1,.5,1)
  rejectc = (.6,.1,.1,1)
badc = (1,0,0,1)
tempc = (.1,1,.1,1)
fixedc = (.9,.9,.9)
count = 0
success = 0
skip = 1
i = 1

colors = [defaultc]*(N+1)
if 'pair' in sys.argv:
  colors[0] = fixedc

def animate(p):
  global count, circles, skip, colors, i, success, plotmoves, contact_counts
  for moves in xrange(plotmoves):
    if (i >= N):
      if 'pair' in sys.argv:
        i = 1
      else:
        i = 0
    keep = True
    temp = move(circles[i])
    if temp[0] - r < 0 or temp[0] + r > lenx:
      keep = False
    else:
      for j in xrange(N):
        if j != i and touch(temp, circles[j]):
          keep = False
          break
    count += 1
    if keep:
      success += 1
    colors[N] = highlightc
    circles[N] = circles[i]

    binnum = round(circles[i][0]/dx) % len(histogram)
    histogram_doubled[2*binnum] += 1 # add another count to histogram
    histogram_doubled[2*binnum+1] += 1 # add another count to histogram
    if keep:
      pairdensity[floor(temp[1]/pairdx), floor(temp[0]/pairdx)] += 1
    else:
      pairdensity[floor(circles[i][1]/pairdx), floor(circles[i][0]/pairdx)] += 1
    if 'density' in sys.argv:
      line.set_ydata(N*area*histogram_doubled/sum(histogram_doubled)/dx/leny)
    oldpos = array([circles[i,0], circles[i,1]])
    if keep:
      circles[i] = temp
    if 'gsigma' in sys.argv:
      contact_counts += 1
      for j in xrange(len(histogram)):
        xjp = j*dx + dx
        xjm = j*dx
        if abs(xjp - circles[i][0]) <= 2 and abs(xjm - circles[i][0]) <= 2:
          lineseg = abs(sqrt(4 - (xjp-circles[i][0])**2) -  sqrt(4 - (xjm-circles[i][0])**2))
          nA_histogram[2*j] += lineseg
          nA_histogram[2*j+1] += lineseg
      inc = zeros(N)
      for j in xrange(N):
        for k in xrange(j):
          if almost_touch(circles[k], circles[j]):
            inc[k] = 1
            inc[j] = 1
      for j in xrange(N):
        if inc[j]:
          binnum = int(round(circles[j][0]/dx)) % len(histogram)
          contact_histogram[2*binnum] += 1
          contact_histogram[2*binnum+1] += 1
    i += 1
  i -= 1

  if keep:
    colors[i] = goodc
  if 'nA' in sys.argv:
    line.set_ydata(100*N*nA_histogram/sum(histogram_doubled))

  ax.cla()
  if 'pair' in sys.argv or 'gsigma' in sys.argv or 'density' in sys.argv:
    # The following is a hokey heuristic to accelerate the clock as we go.
    if str(count).count('0') == len(str(count))-1 and count > 9*plotmoves:
      plotmoves = int(count)
  if 'pair' in sys.argv:
    pairplot = pcolormesh(XXX, YYY, -pairdensity/sum(pairdensity), cmap='hot')
    ax.set_xlim(-edge, lenx+edge)
    ax.set_ylim(-edge, leny+edge)
  if plotmoves == 1:
    circleplotnum = N+1
  else:
    circleplotnum = N
    colors[i] = defaultc
  if plotmoves > 1:
    colors = [defaultc]*(N+1)
    if 'pair' in sys.argv:
      colors[0] = fixedc
  else:
    if not keep:
      colors[N] = badc
      circles[N] = temp
      colors[i] = highlightc
  if 'gsigma' in sys.argv:
    #colors = [defaultc]*(N+1)
    for j in xrange(N):
      for k in xrange(j):
        if almost_touch(circles[k], circles[j]):
          colors[k] = contactc
          colors[j] = contactc
  for j in xrange(circleplotnum):
    for shift in [0, -leny, leny]:
      if 'pair' in sys.argv and j != 0:
        circ = Circle(circles[j]-(0, shift), r, color=colors[j], fill=False, linewidth=3)
      else:
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
  if plotmoves == 1:
    for shift in [0, -leny, leny]:
      delta = periodic_diff(oldpos, temp)
      deltax = temp[0] - oldpos[0]
      deltay = temp[1] - oldpos[1]
      if shift == 0:
        alpha = 1
      else:
        alpha = 0.5
      ax.arrow(oldpos[0], oldpos[1]+shift, delta[0],
               delta[1], head_width=arwidth,
               head_length=arlen, linewidth=2, facecolor='gray', alpha = alpha)
    if keep:
      colors[i] = goodc
    else:
      colors[i] = rejectc
    if (i == N-1):
      colors = [defaultc]*(N+1)
      if 'pair' in sys.argv:
        colors[0] = fixedc
  if 'nA' in sys.argv:
    line.set_ydata(nA_histogram/sum(histogram_doubled))
  if 'gsigma' in sys.argv:
    # FIXME: off by a constant factor below
    gsigma = sum(histogram_doubled)*sum(nA_histogram)*contact_histogram/(histogram_doubled+0.00001)/nA_histogram/delta_contact/80000.0/contact_counts/13.0
    line.set_ydata(gsigma)
  if show:
    time.sleep(1.001)
  i += 1

ax.set_aspect('equal')
ax.set_ylim(-2*edge, leny+2*edge)

for p in xrange(500):
  setup()

for i in xrange(N):
  binnum = round(circles[i][0]/dx) % len(histogram)
  histogram_doubled[2*binnum] += 1 # add another count to histogram
  histogram_doubled[2*binnum+1] += 1 # add another count to histogram

def mysavefig(f):
  savefig(f, transparent=True)
  print 'saved', f, 'with plotmoves', plotmoves, 'and count', count

if 'show' in sys.argv:
  anim = animation.FuncAnimation(fig, animate, init_func=init)
  pylab.show()
else:
  count = 0
  success = 0
  #print('Saving animation images.')
  for p in xrange(31):
    animate(p)
    if 'density' in sys.argv:
      mysavefig("anim/mc-density-%03i.pdf" %p)
    elif 'gsigma' in sys.argv:
      mysavefig("anim/mc-gsigma-%03i.pdf" %p)
    elif 'pair' in sys.argv:
      mysavefig("anim/mc-pair-%03i.pdf" %p)
    else:
      mysavefig("anim/mc-slow-%03i.pdf" %p)

