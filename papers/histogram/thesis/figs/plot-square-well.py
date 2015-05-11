#!/usr/bin/env python2

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plot_dims = (5,4)
plot_scale = 0.7

ww = 2.3
r_max = 4
v_max = 3
v_min = -2

top = 0.3
bottom = -1

plot_size = (plot_dims[0]*plot_scale,plot_dims[1]*plot_scale)

fig = plt.figure(figsize=plot_size)
ax = fig.add_axes([0.15, 0.15, 0.83, 0.83])


plt.plot([1,1,ww,ww,r_max],[v_max,bottom,bottom,top,top],'k')
plt.xlim(0,r_max)
plt.ylim(v_min,v_max)

plt.xticks([0,1,ww],['$0$','$\\sigma$','$\\lambda\\sigma$'])
plt.yticks([top,bottom],['$0$','$-\epsilon$'])
plt.xlabel('$r$')
plt.ylabel('$v_{sw}$')

plt.savefig('figs/square-well.pdf')

