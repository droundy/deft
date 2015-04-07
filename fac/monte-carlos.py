#!/usr/bin/python3

import facfile
import string

src = facfile.facfile('.monte-carlos.fac')

monte_carlos = """
   monte-carlo soft-monte-carlo pair-monte-carlo
   triplet-monte-carlo polyhedra-monte-carlo polyhedra-talk
   square-well-monte-carlo
   radial-distribution-monte-carlo
   free-energy-monte-carlo free-energy-monte-carlo-infinite-case
""".split()

utility_files = """
   polyhedra square-well
""".split()

for x in utility_files:
    src.compile('src/Monte-Carlo/%s.cpp' % x)

for x in monte_carlos:
    src.link('src/Monte-Carlo/%s.cpp' % x, x)
