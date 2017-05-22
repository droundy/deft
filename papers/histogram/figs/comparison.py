from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

T1 = np.genfromtxt(fname = "../deft/papers/histogram/data/s002/periodic-ww1.30-ff0.30-N200-toe-golden-movie/000556-lndos.dat")
print T1



def e_lndos(fbase):
    e_lndos = numpy.loadtxt(fbase+"-lndos.dat", ndmin=2, dtype=numpy.float)

    energy = -e_lndos[:,0] # array of energies
    lndos = e_lndos[:,1]
    return energy, lndos

print e_lndos(
