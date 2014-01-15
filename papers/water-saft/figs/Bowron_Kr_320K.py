#!/usr/bin/python

from __future__ import division

import pylab

# Bowron PRL 1997

bottom_left = (125.0, 755.0)
bottom_left_values = (2.0, 0) # angstroms, dimensionless g

top_right = (907.0, 28.0)
top_right_values = (7.0, 2.5) # angstroms, dimensionless g

pixels = pylab.array([[291.0, 750],
                      [302, 710],
                      [313, 633],
                      [323, 532],
                      [335, 441],
                      [345, 369],
                      [357, 319],
                      [368, 270],
                      [378, 221],
                      [390, 184],
                      [400, 158],
                      [421, 155],
                      [444, 189],
                      [466, 246],
                      [488, 309],
                      [510, 372],
                      [532, 429],
                      [551, 479],
                      [568, 512],
                      [581, 537],
                      [602, 551],
                      [636, 561],
                      [656, 560],
                      [681, 542],
                      [712, 509],
                      [740, 482],
                      [771, 459],
                      [795, 447],
                      [826, 439],
                      [868, 439],
                      [907, 441]])

points = pylab.zeros_like(pixels)

for i in xrange(len(pixels)):
    points[i,0] = (pixels[i,0] - bottom_left[0])*(top_right_values[0] - bottom_left_values[0])/(top_right[0]-bottom_left[0]) + bottom_left_values[0]
    points[i,1] = (pixels[i,1] - bottom_left[1])*(top_right_values[1] - bottom_left_values[1])/(top_right[1]-bottom_left[1]) + bottom_left_values[1]
    #print points[i, 0], points[i, 1]

r = points[:, 0]
g = points[:, 1]

#pylab.plot(r, g)

#pylab.show()
