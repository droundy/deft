#!/usr/bin/env python2

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy

plot_dims = (5,4)
plot_scale = 0.7

ww = 2.3
r_max = 4
v_max = 3
v_min = -2

top = 0.3
bottom = -1

plot_size = (plot_dims[0]*plot_scale,plot_dims[1]*plot_scale)

axis_dimensions = [0.15, 0.15, 0.83, 0.83]

def make_sw_plot():

    fig = plt.figure(figsize=plot_size)
    ax = fig.add_axes(axis_dimensions)

    plt.plot([1,1,ww,ww,r_max],[v_max,bottom,bottom,top,top],'k')
    plt.xlim(0,r_max)
    plt.ylim(v_min,v_max)

    plt.xticks([0,1,ww],['$0$','$\\sigma$','$\\lambda\\sigma$'])
    plt.yticks([top,bottom],['$0$','$-\epsilon$'])
    plt.xlabel('$r$')
    plt.ylabel('$v_{sw}$')

make_sw_plot()
plt.savefig('figs/square-well.pdf')

line_width = 4
color = 'r'

make_sw_plot()
plt.plot([1,1],[v_max,top], color, linewidth=line_width)
plt.savefig('figs/square-well-v1.pdf')

make_sw_plot()
plt.plot([1,1,ww,ww],[top,bottom,bottom,top], color, linewidth=line_width)
plt.savefig('figs/square-well-v2.pdf')

make_sw_plot()
plt.plot([ww,r_max],[top,top], color, linewidth=line_width)
plt.savefig('figs/square-well-v3.pdf')
