#!/usr/bin/python

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab, numpy, sys, scipy.ndimage
import os.path
import math


# these are the things to set
colors = ['k', 'b', 'r', 'g']
plots = ['mc', 'this-work', 'fischer', 'gloor']
plots_dft = ['mc', 'this-work', 'fischer', 'gloor']
dx = 0.1
############################

able_to_read_file = True
z0 = 0.05

#Careful there may be a difference between mc and other with these
zmax = 20
rmax = 10


if len(sys.argv) != 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])


def read_walls(ff, z0, fun):
    global able_to_read_file
    if fun == 'mc':
        filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            able_to_read_file = False
            print "File not found: ", filename
            return 0
    else:
        filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            able_to_read_file = False
            print "File not found: ", filename
            return 0
    print 'Using', filename
    return data


def read_walls_path(ff,z0,fun):
  global able_to_read_file
  if fun == 'mc':
    filename = "figs/mc/wallsMC-pair-%1.1f-path.dat" % ff
    try:
      data = numpy.loadtxt(filename)
    except IOError:
      able_to_read_file = False
      print "File not found: ", filename
      return 0
    print 'Using', filename
    return data[:,0:2]
  else:
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
    try:
      data = numpy.loadtxt(filename)
    except IOError:
      able_to_read_file = False
      print "File not found: ", filename
      return 0
    print 'Using', filename
    return data


g2 = [0]*len(plots)
g2_path = [0]*(len(plots_dft))

num = 100

line = [0]*len(plots)
line_path = [0]*len(plots_dft)
fig = plt.figure()
ax = [0]*2
ax[0] = fig.add_subplot(2,2,1)
ax[1] = fig.add_subplot(2,2,2)

def plot():
    global ax, line
    for i in range(len(plots)):
        g2 = read_walls(ff, z0, plots[i])
        if able_to_read_file == False:
            matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
            matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
            savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
            pylab.savefig(savedfilename)
            continue
        rmax = len(g2[:,0])*dx
        r_array = numpy.arange(0, rmax, dx)
        g2_against_wall = g2[:,0]

        zmax = len(g2[0,:])*dx
        z_array = numpy.arange(dx/2, zmax+dx/2, dx)
        g2_z_direction = g2[0,:]

        theta = numpy.linspace(0, numpy.pi/2, num)
        g2_of_theta = numpy.zeros(len(theta))
        for j in range(len(theta)):
            radius_path=2.5
            r = radius_path*numpy.cos(theta[j])
            z = radius_path*numpy.sin(theta[j])
            g2_rz = g2[numpy.floor(r/dx),numpy.floor(z/dx)]
            g2_rZ = g2[numpy.floor(r/dx),numpy.ceil(z/dx)]
            g2_Rz = g2[numpy.ceil(r/dx),numpy.floor(z/dx)]
            g2_RZ = g2[numpy.ceil(r/dx),numpy.ceil(z/dx)]
            alpha_r = (r/dx)-numpy.floor(r/dx)
            alpha_z = (z/dx)-numpy.floor(z/dx)
            g2_Z = alpha_r*g2_RZ + (1.0-alpha_r)*g2_rZ
            g2_z = alpha_r*g2_Rz + (1.0-alpha_r)*g2_rz
            g2_of_theta[j] = alpha_z*g2_Z + (1.0-alpha_z)*g2_z
        # g2_against_wall = g2_against_wall[::-1] #reverse order of array
        # if plots[i] == 'mc':
        #     r_array = r_array - 5.0
        # for k in numpy.arange(22):
        #     g2_against_wall = g2_against_wall[:len(g2_against_wall)-1]
        #     g2_z_direction = g2_z_direction[1:]
        #     r_array = r_array[:len(r_array)-1]
        #     z_array = z_array[1:]
        # theta = theta + max(r_array)
        # z_array = z_array + max(theta)
        pylab.plot(-r_array[r_array >= radius_path], g2_against_wall[r_array >= radius_path], colors[i], label=plots[i])
        pylab.plot(radius_path*theta - radius_path, g2_of_theta, colors[i])
        pylab.plot((numpy.pi/2-2)*radius_path + z_array[z_array>=radius_path], g2_z_direction[z_array>=radius_path], colors[i])
        #x_axis = numpy.concatenate((r_array,theta,z_array))
        #y_axis = numpy.concatenate((g2_against_wall,g2_of_theta,g2_z_direction))
        #pylab.plot(x_axis,y_axis)
        #pylab.ylim([0,2.5])
        #pylab.xlim([1,18])



def plot_dft():
    global ax, line
    for i in range(1,len(plots)):
        print i
        g2_path = read_walls_path(ff, z0, plots[i])
        if able_to_read_file == False:
            matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
            matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
            savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
            pylab.savefig(savedfilename)
            exit(0)
        x_path = g2_path[:,0]
        pair_correlation = g2_path[:,1]
        pylab.plot(x_path,pair_correlation, colors[i]+'--', label=plots[i])
        #ax[1].set_ylim([0,2.5])
        #ax[1].set_xlim([1,18])

fig = plt.figure(figsize=(10,7))
ax = [0]*2

#fig.add_subplot(2,1,1)
plot()
pylab.xlim(-10,10)
#pylab.legend(loc='best')

#fig.add_subplot(2,1,2)
plot_dft()
pylab.legend(loc='best')
savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
fig.savefig(savedfilename)
pylab.show()
