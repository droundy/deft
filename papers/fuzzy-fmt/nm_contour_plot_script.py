#!/usr/bin/python2
#contour plot

import os
import numpy as np
import matplotlib.pyplot as plt
import sys


fvs = nmdata.variables['fv'][:]
gws = nmdata.variables['gw'][:]
fes = nmdata.variables['fe'][0]
nfvs = lens[fvs]
ngws = lens[gws]


contour_levels = np.arrange()
fvs, gws = np.meshgrid(fvs, gws)
ax = fig.add_axes([0, 1, 0.01, 2])
cs = m.contour()
