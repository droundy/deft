#!/usr/bin/python2
from __future__ import division
import matplotlib, sys, os, argparse
from matplotlib.collections import PolyCollection

parser = argparse.ArgumentParser(description='Plot density and stuff')
parser.add_argument('ff', metavar='ff', type=float, help='filling fraction')
parser.add_argument('-s', '--show', action='store_true')
args = parser.parse_args()

if not args.show:
  matplotlib.use('Agg')
from pylab import *

ff = args.ff
if ff == .3: N = 2400
elif ff == .6: N = 4800
elif ff == .73: N = 5832
else: exit(1)

datdir = "../../papers/polyhedra/figs/mc/"
denfile = datdir + "walls-%4.2f-density-cube-%i.dat" %(ff, N)
orderfile = datdir + "walls-%4.2f-order-cube-%i.dat" %(ff, N)
dendata = loadtxt(denfile)
z = dendata[:,0]
density = dendata[:,3]

order_parameters = transpose(loadtxt(orderfile))
dcostheta = 1/len(order_parameters[:,0])
costheta = arange(0, 1, dcostheta)
Z, Costheta = meshgrid(z, costheta)



fig = figure()
ax = fig.add_subplot(221)
den_ax = fig.add_subplot(222)
order_ax = fig.add_subplot(223)
ax5 = fig.add_subplot(224, sharey = order_ax)
ax5.axis('off')
ax.axis('off')

ax.set_title("Cubes, $ff = %4.2f$" %ff)
den_ax.set_title("number density, $n$")
order_ax.set_title("")


ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
den_ax.get_xaxis().set_visible(False)
den_ax.get_yaxis().set_visible(False)
order_ax.get_xaxis().set_visible(False)

ax5.get_xaxis().set_visible(False)
ax5.get_yaxis().set_visible(False)

ax5.set_aspect('equal')

img = matplotlib.image.imread('figs/cube-img-%4.2f.png' %ff)
ax.imshow(img)

order_ax.set_yticks([cos(arctan(sqrt(2))), 1/sqrt(2), 1-dcostheta])
order_ax.set_yticklabels(['   .   ']*3)

den_ax.plot(z, density)
highest = nanmax(order_parameters.flat)

ax5.set_xlim(0, 1)
#little cubes
scale = .03
x = -.165
cube0 = array([[-1, -1], [-1,1], [1,1], [1,-1]])/2*scale
cube0[:,0] += x
cube0[:,1] += 1-dcostheta
cube1 = array([[0,sqrt(2)], [-sqrt(2), 0], [0, -sqrt(2)], [sqrt(2), 0]])/2*scale
cube1[:,0] += x
cube1[:,1] += 1/sqrt(2)

coll = PolyCollection([cube0, cube1], clip_on=False, color=cm.jet(.3))
ax5.add_collection(coll)

order_ax.set_ylim(1/sqrt(2), 1-dcostheta)

# draw how order parameter is gotten
x = .45
arlen = .05
wally = [0,1]
ax5.plot([x,x], wally, '-k', lw=3)
ax5.arrow(x, .8, -arlen*sqrt(2), 0)

cube = array([[0,sqrt(2)], [-sqrt(2), 0], [0, -sqrt(2)], [sqrt(2), 0]])/2*scale
cube[:,0] += x + .05
cube[:,1] += .8
coll = PolyCollection([cube], color=cm.jet(.3), clip_on=False)

ax5.arrow((cube[0,0] + cube[1,0])/2, (cube[0,1] + cube[1,1])/2, -arlen, arlen)

ax5.add_collection(coll)

order_ax.set_title("rotation")
#order_ax.set_ylabel("angle")

orderplot = order_ax.pcolor(z, costheta, order_parameters, cmap=cm.hot_r, vmax=highest/15, vmin=0)
# levels = linspace(0, highest/10, 20)

# orderplot = order_ax.contourf(z, costheta, order_parameters, levels=levels, extend="max")
cbar_ax = fig.add_axes([0.517, 0.05, 0.02, 0.4])
cs = fig.colorbar(orderplot, cax=cbar_ax, ticks=[])
cs.cmap.set_over('k')
#cs.set_clim([0, highest/10])

fig.tight_layout()

savefig("figs/cube-density-%i.png" %(ff*100))

show()
