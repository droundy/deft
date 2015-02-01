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
  new/Minimize new/NewFunctional
""".split()

for x in generic_sources:
    src.compile('src/%s.cpp' % x)


mkdats = """

  papers/fuzzy-fmt/figs/homogeneous
  papers/fuzzy-fmt/figs/radial-wca
  papers/fuzzy-fmt/figs/soft-wall
  papers/fuzzy-fmt/figs/new-soft-wall
  papers/fuzzy-fmt/figs/new-melting
  papers/fuzzy-fmt/figs/walls
  papers/fuzzy-fmt/figs/new-walls
  papers/fuzzy-fmt/figs/soft-sphere

""".split()

for s in mkdats:
    src.compile(s+'.cpp')
    src.link(s+'.cpp', s+'.mkdat',
             objects = {'src/%s.o' % x for x in generic_sources})
