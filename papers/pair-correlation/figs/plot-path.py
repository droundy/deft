#!/usr/bin/python

from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import pylab, numpy, sys, scipy.ndimage
import os.path
import math


# these are the things to set
colors = ['k', 'b', 'r', 'g']
plots = ['mc', 'this-work', 'fischer', 'gross']
plots_dft = ['this-work', 'fischer', 'gross']
dx = 0.1
############################

able_to_read_file = True

z0 = 0.05
ff = 0.3
dx =.1

#Careful there may be a difference between mc and other with these
zmax = 20
rmax = 10


if len(sys.argv) != 2:
    print("Usage:  " + sys.argv[0] + " ff")
    exit(1)
ff = float(sys.argv[1])


def read_walls(ff, z0, fun):
    if fun == 'mc':
        filename = "figs/mc/wallsMC-pair-%1.1f-%1.2f.dat" % (ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            global able_to_read_file
            able_to_read_file = False
            return 0
    else:
        filename = "figs/walls/wallsWB-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
        try:
            data = numpy.loadtxt(filename)
        except IOError:
            global able_to_read_file
            able_to_read_file = False
            return 0
    print 'Using', filename
    return data


def read_walls_path(ff,z0,fun):
    filename = "figs/walls/wallsWB-path-%s-pair-%1.2f-%1.2f.dat" %(fun, ff, z0)
    try:
        data = numpy.loadtxt(filename)
    except IOError:
        global able_to_read_file
        able_to_read_file = False
        return 0
    print 'Using', filename
    return data


g2 = [0]*len(plots)
g2_path = [0]*(len(plots_dft))

num = 100

line = [0]*4
line_path = [0]*3
fig = plt.figure()
ax = [0]*2
ax[0] = fig.add_subplot(2,2,1)
ax[1] = fig.add_subplot(2,2,2)

def plot():
    global g2, ax, line
    g2[0] = read_walls(ff, z0, 'mc')
    for i in numpy.arange(len(g2)):
        g2[i] = read_walls(ff, z0, plots[i])
        if able_to_read_file == False:
            matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
            matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
            savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
            pylab.savefig(savedfilename)
            exit(0)
        rmax = len(g2[i][:,0])*dx
        r_array = numpy.arange(0, rmax, dx)
        g2_r_array = g2[i][:,0]

        zmax = len(g2[i][0,:])*dx
        z_array = numpy.arange(0, zmax, dx)
        g2_z_array = g2[i][0,:]

        theta = numpy.linspace(0, numpy.pi/2, num)
        g2_of_theta = numpy.zeros(len(theta))
        for j in numpy.arange(len(theta)):
            r = 2.2*numpy.sin(theta[j])
            z = 2.2*numpy.cos(theta[j])
            g2_rz = g2[i][math.floor(r/dx),math.floor(z/dx)]
            g2_rZ = g2[i][math.floor(r/dx),math.ceil(z/dx)]
            g2_Rz = g2[i][math.ceil(r/dx),math.floor(z/dx)]
            g2_RZ = g2[i][math.ceil(r/dx),math.ceil(z/dx)]
            alpha_r = (r/dx)-math.floor(r/dx)
            alpha_z = (z/dx)-math.floor(z/dx)
            g2_Z = alpha_r*g2_RZ + (1.0-alpha_r)*g2_rZ
            g2_z = alpha_r*g2_Rz + (1.0-alpha_r)*g2_rz
            g2_of_theta[j] = alpha_z*g2_Z + (1.0-alpha_z)*g2_z
        g2_r_array = g2_r_array[::-1]
        if i == 0:
            r_array = r_array - 5.0
        for k in numpy.arange(22):
            g2_r_array = g2_r_array[:len(g2_r_array)-1]
            g2_z_array = g2_z_array[1:]
            r_array = r_array[:len(r_array)-1]
            z_array = z_array[1:]
        theta = theta + max(r_array)
        z_array = z_array + max(theta)
        x_axis = numpy.concatenate((r_array,theta,z_array))
        y_axis = numpy.concatenate((g2_r_array,g2_of_theta,g2_z_array))
        line[i], = ax[0].plot(x_axis,y_axis)
        ax[0].set_ylim([0,2.5])
        ax[0].set_xlim([1,18])



def plot_dft():
    global g2_path, ax, line
    for i in numpy.arange(len(g2_path)):
        print i
        g2_path[i] = read_walls_path(ff, z0, plots_dft[i])
        if able_to_read_file == False:
            matplotlib.pyplot.plot(numpy.arange(0,10,1), [0]*10, 'k')
            matplotlib.pyplot.suptitle('!!!!WARNING!!!!! There is data missing from this plot!', fontsize=25)
            savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
            pylab.savefig(savedfilename)
            exit(0)
        x_axis = g2_path[i][:,0]
        y_axis = g2_path[i][:,1]
        line_path[i], = ax[1].plot(x_axis,y_axis)
        ax[1].set_ylim([0,2.5])
        ax[1].set_xlim([1,18])

plot()
plot_dft()
ax[0].legend(line,plots,bbox_to_anchor=(0., 1.3, 1.5, .102))
ax[1].legend(line_path,plots_dft,bbox_to_anchor=(0., 1.3, 1.5, .102))
savedfilename = "figs/pair-correlation-path-" + str(int(ff*10)) + ".pdf"
fig.savefig(savedfilename)
fig.show()

