#!/usr/bin/python

from pylab import *
import os
# Trims the path data to be used for git storage

for eta in [0.1, 0.2, 0.3, 0.4]:
  # pair
  filename = "figs/mc/wallsMC-pair-%03.1f-path.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[800:1900,0:4]
    savetxt("figs/mc/wallsMC-pair-%03.1f-path-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # pair test-particle
  filename = "figs/wallsWB-with-sphere-path-%1.2f.dat" % eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[800:1900,0:4]
    savetxt("figs/wallsWB-with-sphere-path-%1.2f-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # pair 2d:
  filename = "figs/mc/wallsMC-pair-%03.1f-0.01.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[:50, :70]
    savetxt("figs/mc/wallsMC-pair-%03.1f-0.01-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # pair test-particle 2d:
  filename = "figs/wallsWB-with-sphere-%1.2f.dat" % eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[:50, :70]
    savetxt("figs/wallsWB-with-sphere-%1.2f-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # triplet @ contact:
  filename = "figs/mc/triplet/tripletMC-%03.1f-path.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[700:2400, 0:4]
    savetxt("figs/mc/triplet/tripletMC-%3.1f-path-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # triplet @ contact 2d:
  filename = "figs/mc/triplet/tripletMC-%03.1f-02.10.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[:50, :90]
    savetxt("figs/mc/triplet/tripletMC-%03.1f-02.10-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # triplet inbetween
  filename = "figs/mc/triplet/tripletMC-%03.1f-path2.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[400:2500, 0:4]
    savetxt("figs/mc/triplet/tripletMC-%3.1f-path2-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)

  # triplet inbetween 2d:
  filename = "figs/mc/triplet/tripletMC-%03.1f-04.20.dat" %eta
  if (os.path.isfile(filename) == False):
    print("File does not exist: %s" %filename)
  else:
    data = loadtxt(filename)
    trimmed_data = data[:50, :91]
    savetxt("figs/mc/triplet/tripletMC-%03.1f-04.20-trimmed.dat" %eta, trimmed_data, fmt='%.3f')
    print("Trimmed: %s." %filename)
