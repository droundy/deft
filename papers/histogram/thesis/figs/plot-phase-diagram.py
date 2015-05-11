#!/usr/bin/env python2

import matplotlib
matplotlib.use('Agg')
import numpy
import matplotlib.pyplot as plt

import matplotlib.patches as patches

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

plot_dims = (5,4)
plot_scale = 0.7

plot_size = (plot_dims[0]*plot_scale,plot_dims[1]*plot_scale)
fig = plt.figure(figsize=plot_size)
ax = fig.add_axes([0.15, 0.16, 0.8, 0.8])

def draw(start,end,coeffs,color):
    xs = numpy.linspace(0,1)
    ys = coeffs[0]*xs*xs + coeffs[1]*xs
    plt.plot(xs*(end[0]-start[0])+start[0],ys*(end[1]-start[1])+start[1],color)

sg_start = (0,0)
sg_end = (0.3,0.3)
sg_coeffs = (1,0)
draw(sg_start,sg_end,sg_coeffs,'m')

lg_end = (0.6,0.5)
lg_coeffs = (1,0)
draw(sg_end,lg_end,lg_coeffs,'g')

sl_end = (sg_end[0]+0.1,1)
sl_coeffs = (2,2)
draw(sg_end,sl_end,sl_coeffs,'b')

plt.plot(sg_end[0],sg_end[1],'ko')
plt.plot(lg_end[0],lg_end[1],'ko')

def add_text(x,y,text):
    ax.text(x,y,text,ha='center',va='center',
            fontsize=12, color='black',
            transform=ax.transAxes)

add_text(0.15,0.6,'solid')
add_text(0.45,0.5,'liquid')
add_text(0.7,0.2,'gas')
add_text(0.75,0.75,'supercritical\nfluid')

cp_offset = (0.13,0)
add_text(lg_end[0]+cp_offset[0],lg_end[1]+cp_offset[1],
         text='critical\npoint')
tp_offset = (0.1,-0.1)
add_text(sg_end[0]+tp_offset[0],sg_end[1]+tp_offset[1],
         text='triple\npoint')

plt.xlim(0,1)
plt.ylim(0,1)
plt.xticks([0,sg_end[0],lg_end[0]],['$0$','$T_{tp}$','$T_{cr}$'])
plt.yticks([0,sg_end[1],lg_end[1]],['$0$','$P_{tp}$','$P_{cr}$'])
plt.xlabel('Temperature')
plt.ylabel('Pressure')
plt.savefig('figs/phase-diagram.pdf')

