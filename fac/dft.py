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
  papers/fuzzy-fmt/figs/radial-lj
  papers/fuzzy-fmt/figs/new-radial-lj
  papers/fuzzy-fmt/figs/new-radial-wca
  papers/fuzzy-fmt/figs/soft-wall
  papers/fuzzy-fmt/figs/new-soft-wall
  papers/fuzzy-fmt/figs/new-melting
  papers/fuzzy-fmt/figs/walls
  papers/fuzzy-fmt/figs/new-walls
  papers/fuzzy-fmt/figs/new-bh-walls
  papers/fuzzy-fmt/figs/soft-sphere

  papers/square-well-fluid/figs/new-radial-sw
  papers/square-well-fluid/figs/homogeneous
  papers/square-well-fluid/figs/coexistence

  papers/water-saft/figs/surface-tension
  papers/water-saft/figs/equation-of-state
  papers/water-saft/figs/four-rods-in-water
  papers/water-saft/figs/hughes-lj-atom
  papers/water-saft/figs/hughes-lj-atom-hs-density
  papers/water-saft/figs/hughes-single-rod
  papers/water-saft/figs/hughes-sphere
  papers/water-saft/figs/lj-atom
  papers/water-saft/figs/pressure-with-isotherms
  papers/water-saft/figs/rods-in-water
  papers/water-saft/figs/single-rod
  papers/water-saft/figs/sphere

  papers/pair-correlation/figs/walls

""".split()

for s in mkdats:
    # The extra_requirements below is overkill, since it is only
    # needed in pair-correlation.  But it keeps things easy by just
    # adding this as a universal requirement, and costs us little in
    # efficiency.
    src.link(s+'.cpp', s+'.mkdat',
             objects = {'src/%s.o' % x for x in generic_sources},
             extra_requirements = {'papers/pair-correlation/figs/ghs-analytics.h',
                                   'papers/pair-correlation/figs/short-range-ghs-analytics.h'})
