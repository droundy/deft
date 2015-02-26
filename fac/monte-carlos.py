#!/usr/bin/python3

import facfile
import string

src = facfile.facfile('.monte-carlos.fac')

monte_carlos = """
   monte-carlo soft-monte-carlo pair-monte-carlo
   triplet-monte-carlo polyhedra-monte-carlo polyhedra-talk
   square-well-monte-carlo
   radial-distribution-monte-carlo free-energy-monte-carlo
""".split()

utility_files = """
   polyhedra square-well
""".split()

for x in utility_files:
    src.compile('src/Monte-Carlo/%s.cpp' % x)

for x in monte_carlos:
    src.link('src/Monte-Carlo/%s.cpp' % x, x)


for method in ["nw","bubble_suppression","robustly_optimistic","gaussian",
               "wang_landau","optimized_ensemble", 'tmmc', 'kT 1', 'kT 2', 'kT 0.5', 'kT 0.4']:
    for ww in [1.3, 1.5]:
        for ff in [0.3, 0.8]:
            for N in range(5,31):
                outputs = ["papers/histogram/data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
                           % (ww, ff, N, method.replace(' ',''), postfix)
                           for postfix in ['g', 'E', 'lnw', 'transitions']]
                src.rule('./square-well-monte-carlo --%s --N %d --ff %g --ww %g --iterations 1000000'
                         % (method, N, ff, ww),
                         ['square-well-monte-carlo'], outputs)
