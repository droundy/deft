#!/usr/bin/python3

import facfile
import string

src = facfile.facfile('.dft.fac')

generic_sources = """
  lattice utilities Faddeeva
  GridDescription Grid ReciprocalGrid
  IdealGas ChemicalPotential
  HardSpheres ExternalPotential
  Functional ContactDensity
  Gaussian Pow WaterSaftFast WaterSaft_by_handFast
  EffectivePotentialToDensity
  equation-of-state water-constants
  compute-surface-tension
  Minimizer Downhill
  Precision ConjugateGradient
  QuadraticLineMinimizer SteepestDescent
  vector3d
""".split()

for x in generic_sources:
    src.compile('src/%s.cpp' % x)
