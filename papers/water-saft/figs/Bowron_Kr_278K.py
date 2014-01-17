#!/usr/bin/python

from __future__ import division

import pylab

# Bowron PRL 1998

bottom_left = (51.0, 1026.0)
bottom_left_values = (2.0, 0.0) # angstroms, dimensionless g

top_right = (1228.0, 41.0)
top_right_values = (8.0, 3.0) # angstroms, dimensionless g

pixels = pylab.array([[184.0, 1024.0],
                      [226, 1022],
                      [247, 1013],
                      [262, 985],
                      [274, 944],
                      [283, 896],
                      [292, 843],
                      [300, 793],
                      [312, 711],
                      [320, 650],
                      [329, 586],
                      [337, 528],
                      [346, 464],
                      [356, 403],
                      [364, 360],
                      [378, 306],
                      [387, 280],
                      [400, 265],
                      [412, 268],
                      [425, 287],
                      [445, 340],
                      [463, 402],
                      [481, 470],
                      [501, 542],
                      [516, 590],
                      [534, 643],
                      [566, 717],
                      [598, 772],
                      [629, 809],
                      [660, 822],
                      [706, 812],
                      [744, 791],
                      [795, 737],
                      [861, 687],
                      [921, 666],
                      [982, 666],
                      [1044, 670],
                      [1112, 685],
                      [1169, 692],
                      [1226, 697]])

points = pylab.zeros_like(pixels)

for i in xrange(len(pixels)):
    points[i,0] = (pixels[i,0] - bottom_left[0])*(top_right_values[0] - bottom_left_values[0])/(top_right[0]-bottom_left[0]) + bottom_left_values[0]
    points[i,1] = (pixels[i,1] - bottom_left[1])*(top_right_values[1] - bottom_left_values[1])/(top_right[1]-bottom_left[1]) + bottom_left_values[1]
    #print points[i, 0], points[i, 1]

r = points[:, 0]
g = points[:, 1]

#pylab.plot(r, g)

#pylab.show()
