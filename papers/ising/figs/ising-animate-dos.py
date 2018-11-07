#!/usr/bin/env python2

from __future__ import division, print_function
import sys, os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# Specify the filepath.
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, dir_path+'/../results/')

# change this to use colors rather than line at some point
# so that each method gets its own color?
colors = ["cyan","red","green"]

animate_what = sys.argv[1]           # DOS or HIST
normalize_dos = bool(sys.argv[2])    # True or False

# Import the data.
methods = sys.argv[3:]

if animate_what == 'DOS':
  datas = []
  max_S_ever = 0
  min_S_ever = 1e100
  max_Sdiff_ever = 0
  for m in methods:
      d = __import__('%s' % (m))
      d.filename = m
      if d.lndos.max() > max_S_ever:
          max_S_ever = d.lndos.max()
      if d.lndos[d.lndos != 0].min() < min_S_ever:
          min_S_ever = d.lndos[d.lndos != 0].min()
      Sdiff = d.lndos[-1].max() - d.lndos[-1][d.lndos[-1] != 0].min()
      if Sdiff > max_Sdiff_ever:
          print('max_Sdiff_ever', Sdiff, 'from', m)
          max_Sdiff_ever = Sdiff
      datas.append(d) # 'ising1-lnw'
  
  fig = plt.figure()
  ax = plt.axes()
  lines = []
  max_frames = 0
  max_t = []
  for i, d in enumerate(datas):
      line, = plt.plot([], [], color=colors[i], lw=2, label=d.filename) # length 1 tuple
      lines.append(line)
      if len(d.t) > max_frames:
          max_frames = len(d.t)
          max_t = d.t
  plt.legend()
  
  title_text = ax.text(0.05, 0.9, 'hello', transform=ax.transAxes)
  
  def init():
      ax.set_xlim(np.min(datas[0].E),np.max(datas[0].E))
      if normalize_dos:
          ax.set_ylim(-max_Sdiff_ever, 0)
      else:
          ax.set_ylim(min_S_ever,max_S_ever)
      plt.xlabel('E')
      plt.ylabel('$S(E)/k_B$')
      for i in range(len(methods)):
         lines[i].set_data([], [])
         #lines[i].set_color(colors[i])
      return line, title_text
  
  def animate(frame):
    title_text.set_text('t = $10^{%.1f}$' % (np.log10(max_t[frame])))
    to_redraw = [title_text]
    for i in range(len(methods)):
      # ensure that both shortest and longest array remain for duration of animation.
      if frame < len(datas[i].t):
          if normalize_dos:
              lines[i].set_data(datas[i].E, datas[i].lndos[frame]-datas[i].lndos[frame].max())
          else:
              lines[i].set_data(datas[i].E, datas[i].lndos[frame])
      to_redraw.append(lines[i])
    return to_redraw
  
  # movie parameters.
  duration = 10.         # number of seconds the saved movie will last.
  frame_interval = 200. # interval spacing in displayed movie.
  
  # blit=True means only re-draw the parts that have changed.
  anim = animation.FuncAnimation(fig, animate, init_func=init,
                                 frames=max_frames, interval=frame_interval, blit=True)
  
  #anim.save('papers/ising/movies/ising_lndos.mp4', fps=len(datas[0].t)/duration, extra_args=['-vcodec', 'libx264'] )
  plt.show()

if animate_what == 'HIST':
  datas = []
  for m in methods:
      d = __import__('%s' % (m))
      d.filename = m
      datas.append(d) # 'ising1-lnw'
  
  fig = plt.figure()
  ax = plt.axes()
  lines = []
  max_frames = 0
  max_t = []
  for i, d in enumerate(datas):
      line, = plt.plot([], [], color=colors[i], lw=2, label=d.filename) # length 1 tuple
      lines.append(line)
      if len(d.t) > max_frames:
          max_frames = len(d.t)
          max_t = d.t
  plt.legend()
  
  title_text = ax.text(0.05, 0.9, 'hello', transform=ax.transAxes)
  
  def init():
      ax.set_xlim(np.min(datas[0].E),np.max(datas[0].E))
      ax.set_ylim(0, d.histogram[-1].max())
      plt.xlabel('E')
      plt.ylabel('Histogram Counts')
      for i in range(len(methods)):
         lines[i].set_data([], [])
      return line, title_text
  
  def animate(frame):
    title_text.set_text('t = $10^{%.1f}$' % (np.log10(max_t[frame])))
    to_redraw = [title_text]
    for i in range(len(methods)):
      # ensure that both shortest and longest array remain for duration of animation.
      if frame < len(datas[i].t):
              lines[i].set_data(datas[i].E, datas[i].histogram[frame])
      to_redraw.append(lines[i])
    return to_redraw
  
  # movie parameters.
  duration = 10.         # number of seconds the saved movie will last.
  frame_interval = 100. # interval spacing in displayed movie.
  
  # blit=True means only re-draw the parts that have changed.
  anim = animation.FuncAnimation(fig, animate, init_func=init,
                                 frames=max_frames, interval=frame_interval, blit=True)
  
  #anim.save('papers/ising/movies/ising_lndos.mp4', fps=len(datas[0].t)/duration, extra_args=['-vcodec', 'libx264'] )
  plt.show()
