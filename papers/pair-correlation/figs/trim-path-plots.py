#!/usr/bin/python

from pylab import *
import os
# Trims the path data to be used for git storage

for eta in [0.1, 0.2, 0.3, 0.4, 0.5]:
  filename = "figs/mc/wallsMC-pair-%02.1f-path.dat" %eta
  if (os.path.isfile(filename) == False):
    print "File does not exist: ", filename
    continue
  data = loadtxt(filename)
  trimmed_data = data[800:1900,0:4]
  savetxt("figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" %eta, trimmed_data, fmt='%g')

